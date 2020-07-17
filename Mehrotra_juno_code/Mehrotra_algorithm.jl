using LinearAlgebra
using LaTeXStrings
using DataFrames
using DataStructures
using QuadraticModels
using Printf
using SparseArrays
using BenchmarkTools
using NLPModels
using LinearOperators
using QPSReader
using SolverTools
using SolverBenchmark
using Plots
using LDLFactorizations


function starting_points(Qrows, Qcols, Qvals, Arows, Acols, Avals, b, c,
                         lvar, uvar, ilow, iupp, irng, J_augm, n_rows, n_cols, Δ_xλ)

    T = eltype(Avals)

    #J_fact = ldlt(Symmetric(J_augm-Diagonal(tmp_diag), :U))
    J_fact = ldl(Symmetric(J_augm, :U))
    J_P = J_fact.P
    Δ_xλ[n_cols+1: end] = b
    Δ_xλ = ldiv!(J_fact, Δ_xλ)
    #init_xλ2 = J_fact \ [c ; zeros(n_rows)]
    x0 = Δ_xλ[1:n_cols]
    λ0 = Δ_xλ[n_cols+1:end]
    s0_l, s0_u = zeros(n_cols), zeros(n_cols)
    Qx, ATλ = zeros(n_cols), zeros(n_cols)
    Qx = mul_Qx_COO!(Qx, Qrows, Qcols, Qvals, x0)
    ATλ = mul_ATλ_COO!(ATλ, Arows, Acols, Avals, λ0)
    dual_val = Qx - ATλ + c
    s0_l[ilow] = @views dual_val[ilow]
    s0_u[iupp] = @views -dual_val[iupp]

    x0_m_lvar = @views x0[ilow] - lvar[ilow]
    uvar_m_x0 = @views uvar[iupp] - x0[iupp]
    if length(ilow) == 0
        δx_l1, δs_l1 = zero(T), zero(T)
    else
        δx_l1 = max(-T(1.5)*minimum(x0_m_lvar), T(1.e0))
        δs_l1 = @views max(-T(1.5)*minimum(s0_l[ilow]), T(1.e-4))
    end

    if length(iupp) == 0
        δx_u1, δs_u1 = zero(T), zero(T)
    else
        δx_u1 = max(-T(1.5)*minimum(uvar_m_x0), T(1.e0))
        δs_u1 = @views max(-T(1.5)*minimum(s0_u[iupp]), T(1.e-4))
    end

    x0_m_lvar .+= δx_l1
    uvar_m_x0 .+= δx_u1
    s0_l1 = @views s0_l[ilow] .+ δs_l1
    s0_u1 = @views s0_u[iupp] .+ δs_u1
    xs_l1, xs_u1 = s0_l1' * x0_m_lvar, s0_u1' * uvar_m_x0
    if length(ilow) == 0
        δx_l2, δs_l2 = zero(T), zero(T)
    else
        δx_l2 = δx_l1 + xs_l1 / sum(s0_l1) / 2
        δs_l2 = @views δs_l1 + xs_l1 / sum(x0_m_lvar) / 2
    end
    if length(iupp) == 0
        δx_u2, δs_u2 = zero(T), zero(T)
    else
        δx_u2 = δx_u1 + xs_u1 / sum(s0_u1) / 2
        δs_u2 = @views δs_u1 + xs_u1 / sum(uvar_m_x0) / 2
    end
    δx = max(δx_l2, δx_u2)
    δs = max(δs_l2, δs_u2)
    x0[ilow] .+= δx
    x0[iupp] .-= δx
    s0_l[ilow] = @views s0_l[ilow] .+ δs
    s0_u[iupp] = @views s0_u[iupp] .+ δs

    @inbounds @simd for i in irng
        if (lvar[i] < x0[i] < uvar[i]) == false
            x0[i] = (lvar[i] + uvar[i]) / 2
        end
    end

    @assert all(x0 .> lvar) && all(x0 .< uvar)
    @assert @views all(s0_l[ilow] .> zero(T)) && all(s0_u[iupp] .> zero(T))

    return x0, λ0, s0_l, s0_u, J_P, Qx, ATλ, x0_m_lvar, uvar_m_x0, Δ_xλ
end



function compute_α_dual(v, dir_v)
    n = length(v)
    T = eltype(v)
    if n == 0
        return one(T)
    end
    α = one(T)
    @inbounds @simd for i=1:n
        if dir_v[i] < zero(T)
            α_new = -v[i] * T(0.999) / dir_v[i]
            if α_new < α
                α = α_new
            end
        end
    end
    return α
end



function compute_α_primal(v, dir_v, lvar, uvar)
    n = length(v)
    T = eltype(v)
    α_l, α_u = one(T), one(T)
    @inbounds @simd for i=1:n
        if dir_v[i] > zero(T)
            α_u_new = (uvar[i] - v[i]) * T(0.999) / dir_v[i]
            if α_u_new < α_u
                α_u = α_u_new
            end
        elseif dir_v[i] < zero(T)
            α_l_new = (lvar[i] - v[i]) * T(0.999) / dir_v[i]
            if α_l_new < α_l
                α_l = α_l_new
            end
        end
    end
    return min(α_l, α_u)
end



function compute_μ(x_m_lvar, uvar_m_x, s_l, s_u, nb_low, nb_upp)
    return (s_l' * x_m_lvar + s_u' * uvar_m_x) / (nb_low + nb_upp)
end


function solve_augmented_system_aff!(J_fact, Δ_aff, Δ_xλ, rc, rb, x_m_lvar, uvar_m_x,
                                     s_l, s_u, ilow, iupp,  n_cols, n_rows, n_low)

    Δ_xλ[1:n_cols] .= .-rc
    Δ_xλ[n_cols+1:end] .= .-rb
    Δ_xλ[ilow] += @views s_l[ilow]
    Δ_xλ[iupp] -= @views s_u[iupp]

    Δ_xλ = ldiv!(J_fact, Δ_xλ)

    Δ_aff[1:n_cols+n_rows] = Δ_xλ
    Δ_aff[n_cols+n_rows+1:n_cols+n_rows+n_low] = @views -s_l[ilow] - s_l[ilow].*Δ_xλ[1:n_cols][ilow]./x_m_lvar
    Δ_aff[n_cols+n_rows+n_low+1:end] = @views -s_u[iupp] + s_u[iupp].*Δ_xλ[1:n_cols][iupp]./uvar_m_x
    return Δ_aff
end

function solve_augmented_system_cc!(J_fact, Δ_cc, Δ_xλ ,Δ_aff, σ, μ, x_m_lvar, uvar_m_x,
                                    rxs_l, rxs_u, s_l, s_u, ilow, iupp, n_cols, n_rows, n_low)


    rxs_l .= @views (-σ*μ .+ Δ_aff[1:n_cols][ilow].*Δ_aff[n_rows+n_cols+1: n_rows+n_cols+n_low])
    rxs_u .= @views σ*μ .+ Δ_aff[1:n_cols][iupp].*Δ_aff[n_rows+n_cols+n_low+1: end]

    Δ_xλ .= zero(eltype(Δ_xλ))
    Δ_xλ[ilow] .+= rxs_l./x_m_lvar
    Δ_xλ[iupp] .+= rxs_u./uvar_m_x

    Δ_xλ = ldiv!(J_fact, Δ_xλ)

    Δ_cc[1:n_cols+n_rows] = Δ_xλ
    Δ_cc[n_cols+n_rows+1:n_cols+n_rows+n_low] .= @views .-(rxs_l.+s_l[ilow].*Δ_xλ[1:n_cols][ilow])./x_m_lvar
    Δ_cc[n_cols+n_rows+n_low+1:end] .= @views (rxs_u.+s_u[iupp].*Δ_xλ[1:n_cols][iupp])./uvar_m_x
    return Δ_cc
end


function check_frontier!(γf, x, x_m_lvar, uvar_m_x, s_l, s_u, Δ, α_pri, α_dual,
                        x_m_lvar_sup, uvar_m_x_sup, s_l_sup, s_u_sup,
                        lvar, uvar, ilow, iupp, ilow_pad, iupp_pad,
                        n_low, n_upp, n_rows, n_cols)
    T = eltype(x)
    x_m_lvar_sup .= @views x_m_lvar .+ α_pri .* Δ[1:n_cols][ilow]
    uvar_m_x_sup .= @views uvar_m_x .- α_pri .* Δ[1:n_cols][iupp]
    s_l_sup .= @views s_l[ilow] .+ α_dual .* Δ[n_rows+n_cols+1: n_rows+n_cols+n_low]
    s_u_sup .= @views s_u[iupp] .+ α_dual .* Δ[n_rows+n_cols+n_low+1: end]
    μ_sup = compute_μ(x_m_lvar_sup, uvar_m_x_sup, s_l_sup, s_u_sup, n_low, n_upp)
    f_pri, f_pri_tmp = zero(T), zero(T)
    f_dual, f_dual_tmp = zero(T), zero(T)
    i_l, i_u = 0, 0

    @inbounds @simd for i=1:n_cols
        if ilow_pad[i] != 0 && iupp_pad[i] != 0
            i_l, i_u = ilow_pad[i], iupp_pad[i]
            if (x_m_lvar_sup[i_l] == zero(T) || uvar_m_x_sup[i_u] == zero(T))
                f_pri_tmp = (2*γf*μ_sup -
                            x_m_lvar[i_l]*(s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l]) -
                            uvar_m_x[i_u]*(s_u_sup[i_u]+α_dual*Δ[n_rows+n_cols+n_low+i_u])) /
                                (α_pri * Δ[i]*
                                (s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l]-s_u_sup[i_u]-α_dual*Δ[n_rows+n_cols+n_low+i_u]))
                if f_pri_tmp < f_pri
                    f_pri = f_pri_tmp
                end
            end
            if (s_l_sup[i_l] == zero(T) || uvar_m_x_sup[i_u] == zero(T))
                f_dual_tmp = (2*γf*μ_sup - s_l[i_l]*x_m_lvar_sup[i_l] - s_u[i_u]*uvar_m_x_sup[i_u]) /
                                (α_dual * (Δ[n_rows+n_cols+i_l]*x_m_lvar_sup[i_l] +
                                     Δ[n_rows+n_cols+n_low+i_u]*uvar_m_x_sup[i_u]))
                if f_dual_tmp < f_dual
                    f_dual = f_dual_tmp
                end
            end
        elseif ilow_pad[i] != 0 && iupp_pad[i] == 0
            i_l = ilow_pad[i]
            if x_m_lvar_sup[i_l] == zero(T)
                f_pri_tmp = (γf*μ_sup - x_m_lvar[i_l]*(s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l])) /
                                (α_pri * Δ[i]*(s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l]))
                if f_pri_tmp < f_pri
                    f_pri = f_pri_tmp
                end
            end
            if s_l_sup[i_l] == zero(T)
                f_dual_tmp = (γf*μ_sup - s_l[i_l]*x_m_lvar_sup[i_l] ) /
                                (α_dual * Δ[n_rows+n_cols+i_l]*x_m_lvar_sup[i_l])
                if f_dual_tmp < f_dual
                    f_dual = f_dual_tmp
                end
            end
        elseif ilow_pad[i] == 0 && iupp_pad[i] != 0
            i_u = iupp_pad[i]
            if uvar_m_x_sup[i_u] == zero(T)
                f_pri_tmp = (γf*μ_sup - uvar_m_x[i_u]*(s_u_sup[i_u]+α_dual*Δ[n_rows+n_cols+n_low+i_u])) /
                                (α_pri * Δ[i]*(-s_u_sup[i_u]-α_dual*Δ[n_rows+n_cols+n_low+i_u]))
                if f_pri_tmp < f_pri
                    f_pri = f_pri_tmp
                end
            end
            if s_u[i_u] == zero(T)
                f_dual_tmp = (γf*μ_sup - s_u[i_u]*uvar_m_x_sup[i_u]) /
                            (α_dual * Δ[n_rows+n_cols+n_low+i_u]*uvar_m_x_sup[i_u])
                if f_dual_tmp < f_dual
                    f_dual = f_dual_tmp
                end
            end
        end
    end

    α_pri *= max(one(T)-γf, f_pri)
    α_dual *= max(one(T)-γf, f_dual)

    return α_pri, α_dual
end


function mul_Qx_COO!(Qx, Qrows, Qcols, Qvals, x)
    # right mutiplication for sparse COO symetric matrix M: res=Mv
    Qx .= zero(eltype(Qx))
    @inbounds @simd for i=1:length(Qcols)
        Qx[Qrows[i]] += Qvals[i] * x[Qcols[i]]
        if Qrows[i] != Qcols[i]
            Qx[Qcols[i]] += Qvals[i]*x[Qrows[i]]
        end
    end
    return Qx
end

function mul_ATλ_COO!(ATλ, Arows, Acols, Avals, λ)
    ATλ .= zero(eltype(ATλ))
    @inbounds @simd for i=1:length(Acols)
        ATλ[Acols[i]] += Avals[i] * λ[Arows[i]]
    end
    return ATλ
end

function mul_Ax_COO!(Ax, Arows, Acols, Avals, x)
    Ax .= zero(eltype(Ax))
    @inbounds @simd for i=1:length(Acols)
        Ax[Arows[i]] += Avals[i] * x[Acols[i]]
    end
    return Ax
end



function get_norm_rc!(v, A_i, Avals, n_v, n)
    T = eltype(v)
    v .= zero(T)
    @inbounds @simd for j=1:n
        if abs(Avals[j]) > v[A_i[j]]
            v[A_i[j]] = abs(Avals[j])
        end
#         v[A_i[j]] += Avals[j]^2  #2-norm
    end

    v = sqrt.(v)
    @inbounds @simd for i=1:n_v
        if v[i] == zero(T)
            v[i] = one(T)
        end
    end
    return v
end

function mul_A_D1_D2!(Arows, Acols, Avals, d1, d2, r, c, n_rows, n_cols, n)
    @inbounds @simd for i=1:n
        Avals[i] /= r[Arows[i]] * c[Acols[i]]
    end
    d1 ./= r
    d2 ./= c
    return Arows, Acols, Avals, d1, d2
end

function mul_Q_D!(Qrows, Qcols, Qvals, d, c, n_cols, n)
    @inbounds @simd for i=1:n
        Qvals[i] /= c[Qrows[i]] * c[Qcols[i]]
    end
    d ./= c
    return Qrows, Qcols, Qvals, d
end

function scaling_Ruiz!(Arows, Acols, Avals, Qrows, Qcols, Qvals, c, b, lvar, uvar,
                       n_rows, n_cols, ϵ; max_iter = 100)
    n = length(Arows)
    T = eltype(Avals)
    d1, d2 = ones(T, n_rows), ones(T, n_cols)
    r_k, c_k = zeros(T, n_rows), zeros(T, n_cols)

    r_k = get_norm_rc!(r_k, Arows, Avals, n_rows, n)
    c_k = get_norm_rc!(c_k, Acols, Avals, n_cols, n)
    convergence = maximum(abs.(one(T) .- r_k)) <= ϵ && maximum(abs.(one(T) .- c_k)) <= ϵ
    Arows, Acols, Avals, d1, d2 = mul_A_D1_D2!(Arows, Acols, Avals, d1, d2,
                                               r_k, c_k, n_rows, n_cols, n)
    k = 1
    while !convergence && k < max_iter
        r_k = get_norm_rc!(r_k, Arows, Avals, n_rows, n)
        c_k = get_norm_rc!(c_k, Acols, Avals, n_cols, n)
        convergence = maximum(abs.(one(T) .- r_k)) <= ϵ && maximum(abs.(one(T) .- c_k)) <= ϵ
        Arows, Acols, Avals, d1, d2 = mul_A_D1_D2!(Arows, Acols, Avals, d1, d2,
                                                   r_k, c_k, n_rows, n_cols, n)
        k += 1
    end

    n_Q = length(Qrows)
    @inbounds @simd for i=1:n_Q
        Qvals[i] *= d2[Qrows[i]] * d2[Qcols[i]]
    end
    b .*= d1
    c .*= d2
    lvar ./= d2
    uvar ./= d2

    # scaling Q (symmetric)
    d3 = ones(T, n_cols)
    c_k .= zero(T)
    c_k = get_norm_rc!(c_k, Qcols, Qvals, n_cols, n_Q)
    convergence = maximum(abs.(one(T) .- c_k)) <= ϵ
    Qrows, Qcols, Qvals, d3 = mul_Q_D!(Qrows, Qcols, Qvals, d3, c_k, n_cols, n_Q)
    k = 1
    while !convergence && k < max_iter
        c_k = get_norm_rc!(c_k, Qcols, Qvals, n_cols, n_Q)
        convergence = maximum(abs.(one(T) .- c_k)) <= ϵ
        Qrows, Qcols, Qvals, d3 = mul_Q_D!(Qrows, Qcols, Qvals, d3, c_k, n_cols, n_Q)
        k += 1
    end

    for i=1:n
        Avals[i] *= d3[Acols[i]]
    end
    c .*= d3
    lvar ./= d3
    uvar ./= d3

    return Arows, Acols, Avals, Qrows, Qcols, Qvals, c, b, lvar, uvar, d1, d2, d3
end


function get_diag_sparseCSC(M; tri=:U)
    # get diagonal index of M.nzval
    # we assume all columns of M are non empty, and M triangular (:L or :U)
    @assert tri ==:U || tri == :L
    T = eltype(M)
    n = length(M.rowval)
    diagind = zeros(Int, M.m) # square matrix
    index = M.rowval[1] # 1
    if tri == :U
        for i=1:M.m
            diagind[i] = M.colptr[i+1] - 1
        end
    else
        for i=1:M.m
            diagind[i] = M.colptr[i]
        end
    end
    return diagind
end

function get_diag_sparseCOO(Qrows, Qcols, Qvals, n_cols)
    # get diagonal index of M.nzval
    # we assume all columns of M are non empty, and M triangular (:L or :U)
    T = eltype(Qvals)
    n = length(Qrows)
    diagval = zeros(T, n_cols)
    for i=1:n
        if Qrows[i] == Qcols[i]
            diagval[Qrows[i]] = Qvals[i]
        end
    end

    return diagval
end

function mehrotraPCQuadBounds(QM0; max_iter=200, ϵ_pdd=1e-8, ϵ_rb=1e-6, ϵ_rc=1e-6,
                              tol_Δx=1e-16, ϵ_μ=0., max_time=1080., scaling=true,
                              display=true, check_frontier=true)
    QM = SlackModel(QM0)

    start_time = time()
    elapsed_time = 0.0

    # get variables from QuadraticModel
    lvar, uvar = QM.meta.lvar, QM.meta.uvar
    n_cols = length(lvar)
    Oc = zeros(n_cols)
    ilow, iupp = [QM.meta.ilow; QM.meta.irng], [QM.meta.iupp; QM.meta.irng] # finite bounds index
    n_low, n_upp = length(ilow), length(iupp) # number of finite constraints
    c = grad(QM, Oc)
    A = jac(QM, Oc)
    T = eltype(A)
    Arows, Acols, Avals = findnz(A)
    n_rows, n_cols = size(A)
    @assert QM.meta.lcon == QM.meta.ucon # equality constraint (Ax=b)
    b = QM.meta.lcon
    Q = hess(QM, Oc)  # lower triangular
    Qrows, Qcols, Qvals = findnz(Q)
    c0 = obj(QM, Oc)

    if scaling
        Arows, Acols, Avals, Qrows, Qcols, Qvals,
        c, b, lvar, uvar, d1, d2, d3 = scaling_Ruiz!(Arows, Acols, Avals, Qrows, Qcols, Qvals,
                                                     c, b, lvar, uvar, n_rows, n_cols, T(1e-4))
    end


    # init regularization values
    ρ, δ = T(1e5*sqrt(eps())), T(1e6*sqrt(eps())) # 1e6, 1e-1 ok

    J_augmrows = vcat(Qcols, Acols, n_cols+1:n_cols+n_rows, 1:n_cols)
    J_augmcols = vcat(Qrows, Arows.+n_cols, n_cols+1:n_cols+n_rows, 1:n_cols)
    tmp_diag = -T(1.0e-7).*ones(T, n_cols)
    J_augmvals = vcat(-Qvals, Avals, δ*ones(n_rows), tmp_diag)
    J_augm = sparse(J_augmrows, J_augmcols, J_augmvals)
    diagind_J = get_diag_sparseCSC(J_augm)
    diag_Q = get_diag_sparseCOO(Qrows, Qcols, Qvals, n_cols)

    k = 0
    Δ_aff = zeros(T, n_cols+n_rows+n_low+n_upp)
    Δ_cc = zeros(T, n_cols+n_rows+n_low+n_upp)
    Δ = zeros(T, n_cols+n_rows+n_low+n_upp)
    Δ_xλ = zeros(T, n_cols+n_rows)

    x, λ, s_l, s_u, J_P, Qx, ATλ,
    x_m_lvar, uvar_m_x, Δ_xλ = @views starting_points(Qrows, Qcols, Qvals, Arows, Acols, Avals,
                                                      b, c, lvar, uvar, ilow, iupp, QM.meta.irng,
                                                      J_augm , n_rows, n_cols, Δ_xλ)

    Qx = mul_Qx_COO!(Qx, Qrows, Qcols, Qvals, x)
    ATλ = mul_ATλ_COO!(ATλ, Arows, Acols, Avals, λ)
    Ax = zeros(n_rows)
    Ax = mul_Ax_COO!(Ax, Arows, Acols, Avals, x)
    rb = Ax - b
    rc = -Qx + ATλ + s_l - s_u - c

    x_m_lvar .= @views x[ilow] .- lvar[ilow]
    uvar_m_x .= @views uvar[iupp] .- x[iupp]
    μ = @views compute_μ(x_m_lvar, uvar_m_x, s_l[ilow], s_u[iupp], n_low, n_upp)

    x_m_l_αΔ_aff = zeros(T, n_low) # x-lvar + αΔ_aff
    u_m_x_αΔ_aff = zeros(T, n_upp) # uvar-x + αΔ_aff
    s_l_αΔ_aff = zeros(T, n_low) # s_l + αΔ_aff
    s_u_αΔ_aff = zeros(T, n_upp) # s_l + αΔ_aff
    rxs_l, rxs_u = zeros(T, n_low), zeros(T, n_upp)
    ilow_pad, iupp_pad = zeros(Int, n_cols), zeros(Int, n_cols)
    ilow_pad[ilow] = 1:n_low
    iupp_pad[iupp] = 1:n_upp

    if check_frontier
        x_m_lvar_sup, uvar_m_x_sup = zeros(T, n_low), zeros(T, n_upp)
        s_l_sup, s_u_sup = zeros(T, n_low), zeros(T, n_upp)
    end

    # stopping criterion
    xTQx_2 = x' * Qx / 2
    cTx = c' * x
    pri_obj = xTQx_2 + cTx + c0
    dual_obj = b' * λ - xTQx_2 + view(s_l,ilow)'*view(lvar,ilow) -
                    view(s_u,iupp)'*view(uvar,iupp) +c0
    pdd = abs(pri_obj - dual_obj ) / (one(T) + abs(pri_obj))
    max_rc, max_rb = norm(rc, Inf), norm(rb, Inf)
    optimal = pdd < ϵ_pdd && max_rb < ϵ_rb && max_rc < ϵ_rc

    n_Δx = zero(T)
    small_Δx, small_μ = false, μ < ϵ_μ
    Δt = time() - start_time
    tired = Δt > max_time

    # display
    if display == true
        @info log_header([:k, :pri_obj, :pdd, :max_rb, :max_rc, :n_Δx, :μ],
                         [Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
                        hdr_override=Dict(:k=>"Iter", :pri_obj=>"primal", :pdd=>"pdd",
                                          :max_rb=>"rb cond", :max_rc=>"rc cond",
                                          :n_Δx=>"‖Δx‖", :μ=>"μ"))
        @info log_row([k, pri_obj, pdd, max_rb, max_rc, n_Δx, μ])
    end

    while k<max_iter && !optimal && !tired # && !small_μ && !small_μ

            # Affine scaling direction
        tmp_diag .= -ρ
        tmp_diag[ilow] .-= @views s_l[ilow] ./ x_m_lvar
        tmp_diag[iupp] .-= @views s_u[iupp] ./ uvar_m_x

#         J_augmvals[end-n_cols+1:end] = tmp_diag
#         J_augm = sparse(J_augmrows, J_augmcols, J_augmvals)

        J_augm.nzval[view(diagind_J,1:n_cols)] .= @views tmp_diag .- diag_Q
        J_augm.nzval[view(diagind_J, n_cols+1:n_rows+n_cols)] .= δ


        J_fact = ldl(Symmetric(J_augm, :U), J_P)

        Δ_aff = solve_augmented_system_aff!(J_fact, Δ_aff, Δ_xλ, rc, rb, x_m_lvar, uvar_m_x,
                                            s_l, s_u, ilow, iupp,  n_cols, n_rows, n_low)

        α_aff_pri = @views compute_α_primal(x, Δ_aff[1:n_cols], lvar, uvar)
        α_aff_dual_l = @views compute_α_dual(s_l[ilow], Δ_aff[n_rows+n_cols+1: n_rows+n_cols+n_low])
        α_aff_dual_u = @views compute_α_dual(s_u[iupp], Δ_aff[n_rows+n_cols+n_low+1:end])

        # alpha_aff_dual_final is the min of the 2 alpha_aff_dual
        α_aff_dual_final = min(α_aff_dual_l, α_aff_dual_u)

        x_m_l_αΔ_aff .= @views x_m_lvar .+ α_aff_pri .* Δ_aff[1:n_cols][ilow]
        u_m_x_αΔ_aff .= @views uvar_m_x .- α_aff_pri .* Δ_aff[1:n_cols][iupp]
        s_l_αΔ_aff .= @views s_l[ilow] .+ α_aff_dual_final .* Δ_aff[n_rows+n_cols+1: n_rows+n_cols+n_low]
        s_u_αΔ_aff .= @views s_u[iupp] .+ α_aff_dual_final .*  Δ_aff[n_rows+n_cols+n_low+1: end]

        μ_aff = compute_μ(x_m_l_αΔ_aff, u_m_x_αΔ_aff, s_l_αΔ_aff, s_u_αΔ_aff,
                          n_low, n_upp)

        σ = (μ_aff / μ)^3

        # corrector and centering step

        Δ_cc = solve_augmented_system_cc!(J_fact, Δ_cc, Δ_xλ , Δ_aff, σ, μ,x_m_lvar, uvar_m_x,
                                          rxs_l, rxs_u, s_l, s_u, ilow, iupp, n_cols, n_rows, n_low)


        Δ .= Δ_aff .+ Δ_cc # final direction

        α_pri = @views compute_α_primal(x, Δ[1:n_cols], lvar, uvar)
        α_dual_l = @views compute_α_dual(s_l[ilow], Δ[n_rows+n_cols+1: n_rows+n_cols+n_low])
        α_dual_u = @views compute_α_dual(s_u[iupp], Δ[n_rows+n_cols+n_low+1: end])

        α_dual_final = min(α_dual_l, α_dual_u)

        if check_frontier
            α_pri, α_dual_final = check_frontier!(T(1.0e-3), x, x_m_lvar, uvar_m_x, s_l, s_u,
                                                  Δ, α_pri, α_dual_final,
                                                  x_m_lvar_sup, uvar_m_x_sup, s_l_sup, s_u_sup,
                                                  lvar, uvar, ilow, iupp,
                                                  ilow_pad, iupp_pad, n_low, n_upp, n_rows, n_cols)
        end

        # new parameters
        x .= @views x .+ α_pri .* Δ[1:n_cols]
        λ .= @views λ .+ α_dual_final .* Δ[n_cols+1: n_rows+n_cols]
        s_l[ilow] .= @views s_l[ilow] .+ α_dual_final .* Δ[n_rows+n_cols+1: n_rows+n_cols+n_low]
        s_u[iupp] .= @views s_u[iupp] .+ α_dual_final .* Δ[n_rows+n_cols+n_low+1: end]
        n_Δx = @views α_pri * norm(Δ[1:n_cols])
        x_m_lvar .= @views x[ilow] .- lvar[ilow]
        uvar_m_x .= @views uvar[iupp] .- x[iupp]

        μ = @views compute_μ(x_m_lvar, uvar_m_x, s_l[ilow], s_u[iupp],
                             n_low, n_upp)

        Qx = mul_Qx_COO!(Qx, Qrows, Qcols, Qvals, x)
        xTQx_2 =  x' * Qx / 2
        ATλ = mul_ATλ_COO!(ATλ, Arows, Acols, Avals, λ)
        Ax = mul_Ax_COO!(Ax, Arows, Acols, Avals, x)
        cTx = c' * x
        pri_obj = xTQx_2 + cTx + c0
        dual_obj = b' * λ - xTQx_2 + view(s_l,ilow)'*view(lvar,ilow) -
                    view(s_u,iupp)'*view(uvar,iupp) +c0

        rb .= Ax .- b
        rc .= ATλ .-Qx .+ s_l .- s_u .- c

        # update stopping criterion values:

        pdd = abs(pri_obj - dual_obj ) / (one(T) + abs(pri_obj))
        max_rc, max_rb = norm(rc, Inf), norm(rb, Inf)
        optimal = pdd < ϵ_pdd && max_rb < ϵ_rb && max_rc < ϵ_rc
        small_Δx, small_μ = n_Δx < tol_Δx, μ < ϵ_μ
        k += 1

        if δ >= 100*T(sqrt(eps()))
            δ /= 10
            #J_augmvals[end-n_cols-n_rows+1:end-n_cols] .= δ
        end
        if ρ >= 1*T(sqrt(eps()))
            ρ /= 10
        end

        Δt = time() - start_time
        tired = Δt > max_time

        if display == true
            @info log_row([k, pri_obj, pdd, max_rb, max_rc, n_Δx, μ])
        end

    end

    if k>= max_iter
        status = :max_iter
    elseif tired
        status = :max_time
    elseif optimal
        status = :acceptable
    else
        status = :unknown
    end

    if scaling
        x .*= d2 .* d3
        for i=1:length(Qrows)
            Qvals[i] /= d2[Qrows[i]] * d2[Qcols[i]] * d3[Qrows[i]] * d3[Qcols[i]]
        end
        Qx = mul_Qx_COO!(Qx, Qrows, Qcols, Qvals, x)
        xTQx_2 =  x' * Qx / 2
        for i=1:length(Arows)
            Avals[i] /= d1[Arows[i]] * d2[Acols[i]] * d3[Acols[i]]
        end
        λ .*= d1
        ATλ = mul_ATλ_COO!(ATλ, Arows, Acols, Avals, λ)
        Ax = mul_Ax_COO!(Ax, Arows, Acols, Avals, x)
        b ./= d1
        c ./= d2 .* d3
        cTx = c' * x
        pri_obj = xTQx_2 + cTx + c0
        lvar .*= d2 .* d3
        uvar .*= d2 .* d3
        dual_obj = b' * λ - xTQx_2 + view(s_l,ilow)'*view(lvar,ilow) -
                    view(s_u,iupp)'*view(uvar,iupp) +c0

        s_l ./= d2 .* d3
        s_u ./= d2 .* d3
        rb .= Ax .- b
        rc .= ATλ .-Qx .+ s_l .- s_u .- c
        max_rc, max_rb = norm(rc, Inf), norm(rb, Inf)
    end

    elapsed_time = time() - start_time

    stats = GenericExecutionStats(status, QM, solution = x,
                                  objective = pri_obj ,
                                  dual_feas = max_rc,
                                  primal_feas = max_rb,
                                  multipliers = λ,
                                  multipliers_L = s_l,
                                  multipliers_U = s_u,
                                  iter = k, elapsed_time=elapsed_time)
    return stats
end




function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end
