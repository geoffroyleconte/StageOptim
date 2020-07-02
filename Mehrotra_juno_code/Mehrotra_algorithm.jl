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


function init_x0_lsq(A, b, lvar, uvar)
    x_tilde = A\b
    n = length(x_tilde)
    for i=1:n
        if x_tilde[i] <= lvar[i]
            x_tilde[i] = lvar[i] + 1
        elseif uvar[i] <= x_tilde[i]
            x_tilde[i] = uvar[i] - 1
        end
        if !(lvar[i] < x_tilde[i] < uvar[i])
            x_tilde[i] = (lvar[i] + uvar[i])/2
        end
    end

    return x_tilde
end

function starting_points3(Q, A, b, c, lvar, uvar, ilow, iupp, irng, J_augm, ρ, n_rows, n_cols)

    tmp_diag = vcat(ρ*ones(n_cols),spzeros(n_rows))

    #J_fact = ldlt(Symmetric(J_augm-Diagonal(tmp_diag), :U))
    J_fact = ldl(Symmetric(J_augm-Diagonal(tmp_diag), :U))
    J_P = J_fact.P
    init_xλ = J_fact \ [zeros(n_cols) ; b]
    #init_xλ2 = J_fact \ [c ; zeros(n_rows)]
    x0 = init_xλ[1:n_cols]
    λ0 = init_xλ[n_cols+1:end]
    s0_l, s0_u = zeros(n_cols), zeros(n_cols)
    dual_val = Q*x0 - A'*λ0 + c
    s0_l[ilow] = @views dual_val[ilow]
    s0_u[iupp] = @views -dual_val[iupp]

    x0_l1, x0_u1 = x0[ilow], x0[iupp]
    if length(ilow) == 0
        δx_l1, δs_l1 = 0., 0.
    else
        δx_l1 = @views max(-1.5*minimum(x0_l1 - lvar[ilow]), 1.)
        δs_l1 = @views max(-1.5*minimum(s0_l[ilow]), 1.e-4)
    end

    if length(iupp) == 0
        δx_u1, δs_u1 = 0., 0.
    else
        δx_u1 = @views max(-1.5*minimum((uvar[iupp] - x0_u1)), 1.)
        δs_u1 = @views max(-1.5*minimum(s0_u[iupp]), 1.e-4)
    end

    x0_l1 .+= δx_l1
    x0_u1 .-= δx_u1
    s0_l1 = @views s0_l[ilow] .+ δs_l1
    s0_u1 = @views s0_u[iupp] .+ δs_u1
    xs_l1, xs_u1 = @views s0_l1'*(x0_l1 - lvar[ilow]), s0_u1'*(uvar[iupp] - x0_u1)
    if length(ilow) == 0
        δx_l2, δs_l2 = 0., 0.
    else
        δx_l2 = δx_l1 + 0.5 *  xs_l1 / sum(s0_l1)
        δs_l2 = @views δs_l1 + 0.5 * xs_l1 / sum(x0_l1-lvar[ilow])
    end
    if length(iupp) == 0
        δx_u2, δs_u2 = 0., 0.
    else
        δx_u2 = δx_u1 + 0.5 *  xs_u1 / sum(s0_u1)
        δs_u2 = @views δs_u1 + 0.5 * xs_u1 / sum(uvar[iupp]-x0_u1)
    end
    δx = max(δx_l2, δx_u2)
    δs = max(δs_l2, δs_u2)
    x0[ilow] .+= δx
    x0[iupp] .-= δx
    s0_l[ilow] = @views s0_l[ilow] .+ δs
    s0_u[iupp] = @views s0_u[iupp] .+ δs

    for i in irng
        if (lvar[i] < x0[i] < uvar[i]) == false
            x0[i] = (lvar[i] + uvar[i]) / 2.
        end
    end

    @assert all(x0 .> lvar) && all(x0 .< uvar)
    @assert @views all(s0_l[ilow] .> 0) && all(s0_u[iupp] .> 0)

    return x0, λ0, s0_l, s0_u, J_P
end


function compute_α_dual(v, dir_v)
    n = length(v)
    if n == 0
        return 1.
    end
    α = 1.
    for i=1:n
        if dir_v[i] < 0.
            α_new = -v[i] * 0.995 / dir_v[i]
            if α_new < α
                α = α_new
            end
        end
    end
    return α
end



function compute_α_primal(v, dir_v, lvar, uvar)
    n = length(v)
    α_l, α_u = 1., 1.
    for i=1:n
        if dir_v[i] > 0.
            α_u_new = (uvar[i] - v[i]) * 0.995 / dir_v[i]
            if α_u_new < α_u
                α_u = α_u_new
            end
        elseif dir_v[i] < 0.
            α_l_new = (lvar[i] - v[i]) * 0.995 / dir_v[i]
            if α_l_new < α_l
                α_l = α_l_new
            end
        end
    end
    return min(α_l, α_u)
end



function compute_μ(x_minus_lvar, uvar_minus_x, s_l, s_u, nb_low, nb_upp)
    return (s_l' * x_minus_lvar + s_u' * uvar_minus_x) / (nb_low + nb_upp)
end


function is_in_Neighborhood_inf(gamma, x_l, x_u, s_l, s_u, lvar, uvar)
# check if the current point is in N_inf(gamma)
# true : (xi_l - lvari) * si_l >= gamma mu   and   (uvari - xi_u) * si_u >= gamma mu
mu = Compute_mu(x_l, x_u, s_l, s_u, lvar, uvar)
for i=1:length(x_l)
    if (x_l[i] - lvar[i]) * s_l[i] < gamma*mu
        return false
    end
end
for i=1:length(x_u)
    if (uvar[i] - x_u[i]) * s_u[i] < gamma*mu
        return false
    end
end
return true
end

function solve_augmented_system_aff!(J_fact, Δ_aff, Δ_x_λ, rc, rb, x_minus_lvar, uvar_minus_x,
                                     s_l, s_u, ilow, iupp,  n_cols, n_rows, n_low)

    F_x_λ_aff = [-rc
                 -rb]
    F_x_λ_aff[ilow] += @views s_l[ilow]
    F_x_λ_aff[iupp] -= @views s_u[iupp]

    Δ_x_λ = ldiv!(Δ_x_λ , J_fact, F_x_λ_aff)
#     Δ_aff = @views [Δ_x_λ
#                     -s_l[ilow] - s_l[ilow].*Δ_x_λ[1:n_cols][ilow]./x_minus_lvar
#                     -s_u[iupp] + s_u[iupp].*Δ_x_λ[1:n_cols][iupp]./uvar_minus_x]
    Δ_aff[1:n_cols+n_rows] = Δ_x_λ
    Δ_aff[n_cols+n_rows+1:n_cols+n_rows+n_low] = @views -s_l[ilow] - s_l[ilow].*Δ_x_λ[1:n_cols][ilow]./x_minus_lvar
    Δ_aff[n_cols+n_rows+n_low+1:end] = @views -s_u[iupp] + s_u[iupp].*Δ_x_λ[1:n_cols][iupp]./uvar_minus_x
    return Δ_aff
end

function solve_augmented_system_cc!(J_fact, Δ_cc, Δ_x_λ ,Δ_aff, σ, μ,x_minus_lvar, uvar_minus_x,
                                    s_l, s_u, ilow, iupp, n_cols, n_rows, n_low)


    rxs_l = @views -σ*μ .+ Δ_aff[1:n_cols][ilow].*Δ_aff[n_rows+n_cols+1: n_rows+n_cols+n_low]
    rxs_u = @views σ*μ .+ Δ_aff[1:n_cols][iupp].*Δ_aff[n_rows+n_cols+n_low+1: end]

    F_x_λ_cc = zeros(n_cols+n_rows)
    F_x_λ_cc[ilow] += rxs_l./x_minus_lvar
    F_x_λ_cc[iupp] += rxs_u./uvar_minus_x

#     F_x_λ_cc = [ I_low*inv_X_L*rxs_l + I_op_upp*inv_X_U*rxs_u
#                                    Or                 ]

    Δ_x_λ = ldiv!(Δ_x_λ , J_fact, F_x_λ_cc)

#     Δ_cc = @views [Δ_x_λ
#                    -(rxs_l+s_l[ilow].*Δ_x_λ[1:n_cols][ilow])./x_minus_lvar
#                    (rxs_u+s_u[iupp].*Δ_x_λ[1:n_cols][iupp])./uvar_minus_x ]
    Δ_cc[1:n_cols+n_rows] = Δ_x_λ
    Δ_cc[n_cols+n_rows+1:n_cols+n_rows+n_low] = @views -(rxs_l+s_l[ilow].*Δ_x_λ[1:n_cols][ilow])./x_minus_lvar
    Δ_cc[n_cols+n_rows+n_low+1:end] = @views (rxs_u+s_u[iupp].*Δ_x_λ[1:n_cols][iupp])./uvar_minus_x
    return Δ_cc
end


function check_frontier!(γf, x, x_minus_lvar, uvar_minus_x, s_l, s_u, Δ, α_pri, α_dual,
                        lvar, uvar, ilow, iupp, n_low, n_upp, n_rows, n_cols)

    x_minus_lvar_sup = @views a_plus_αb(x_minus_lvar, α_pri, Δ[1:n_cols][ilow], n_low)
    uvar_minus_x_sup = @views a_plus_αb(uvar_minus_x, -α_pri, Δ[1:n_cols][iupp], n_upp)
    s_l_sup = @views a_plus_αb(s_l[ilow], α_dual, Δ[n_rows+n_cols+1: n_rows+n_cols+n_low], n_low)
    s_u_sup = @views a_plus_αb(s_u[iupp], α_dual,  Δ[n_rows+n_cols+n_low+1: end], n_upp)
    μ_sup = compute_μ(x_minus_lvar_sup, uvar_minus_x_sup, s_l_sup, s_u_sup, n_low, n_upp)
    f_pri, f_pri_tmp = 0., 0.
    f_dual, f_dual_tmp = 0., 0.
    i_l, i_u = 0, 0
    ilow_pad, iupp_pad = zeros(Int, n_cols), zeros(Int, n_cols)
    ilow_pad[ilow] = 1:n_low
    iupp_pad[iupp] = 1:n_upp
    for i=1:n_cols
        if ilow_pad[i] != 0 && iupp_pad[i] != 0
            i_l, i_u = ilow_pad[i], iupp_pad[i]
            if (x_minus_lvar_sup[i_l] == 0. || uvar_minus_x_sup[i_u] == 0.)
                f_pri_tmp = (2.0*γf*μ_sup -
                            x_minus_lvar[i_l]*(s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l]) -
                            uvar_minus_x[i_u]*(s_u_sup[i_u]+α_dual*Δ[n_rows+n_cols+n_low+i_u])) /
                                (α_pri * Δ[i]*
                                (s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l]-s_u_sup[i_u]-α_dual*Δ[n_rows+n_cols+n_low+i_u]))
                if f_pri_tmp < f_pri
                    f_pri = f_pri_tmp
                end
            end
            if (s_l_sup[i_l] == 0. || uvar_minus_x_sup[i_u] == 0.)
                f_dual_tmp = (2.0*γf*μ_sup - s_l[i_l]*x_minus_lvar_sup[i_l] - s_u[i_u]*uvar_minus_x_sup[i_u]) /
                                (α_dual * (Δ[n_rows+n_cols+i_l]*x_minus_lvar_sup[i_l] +
                                     Δ[n_rows+n_cols+n_low+i_u]*uvar_minus_x_sup[i_u]))
                if f_dual_tmp < f_dual
                    f_dual = f_dual_tmp
                end
            end
        elseif ilow_pad[i] != 0 && iupp_pad[i] == 0
            i_l = ilow_pad[i]
            if x_minus_lvar_sup[i_l] == 0.
                f_pri_tmp = (γf*μ_sup - x_minus_lvar[i_l]*(s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l])) /
                                (α_pri * Δ[i]*(s_l_sup[i_l]+α_dual*Δ[n_rows+n_cols+i_l]))
                if f_pri_tmp < f_pri
                    f_pri = f_pri_tmp
                end
            end
            if s_l_sup[i_l] == 0.
                f_dual_tmp = (γf*μ_sup - s_l[i_l]*x_minus_lvar_sup[i_l] ) /
                                (α_dual * Δ[n_rows+n_cols+i_l]*x_minus_lvar_sup[i_l])
                if f_dual_tmp < f_dual
                    f_dual = f_dual_tmp
                end
            end
        elseif ilow_pad[i] == 0 && iupp_pad[i] != 0
            i_u = iupp_pad[i]
            if uvar_minus_x_sup[i_u] == 0.
                f_pri_tmp = (γf*μ_sup - uvar_minus_x[i_u]*(s_u_sup[i_u]+α_dual*Δ[n_rows+n_cols+n_low+i_u])) /
                                (α_pri * Δ[i]*(-s_u_sup[i_u]-α_dual*Δ[n_rows+n_cols+n_low+i_u]))
                if f_pri_tmp < f_pri
                    f_pri = f_pri_tmp
                end
            end
            if s_u[i_u] == 0.
                f_dual_tmp = (γf*μ_sup - s_u[i_u]*uvar_minus_x_sup[i_u]) /
                            (α_dual * Δ[n_rows+n_cols+n_low+i_u]*uvar_minus_x_sup[i_u])
                if f_dual_tmp < f_dual
                    f_dual = f_dual_tmp
                end
            end
        end
    end

    α_pri *= max(1.0-γf, f_pri)
    α_dual *= max(1.0-γf, f_dual)

    return α_pri, α_dual
end



function a_plus_αb(a, α, b, n)
    # a, b same length, α float
    result = zeros(n)
    for i=1:n
        result[i] = a[i] + α * b[i]
    end
    return result
end

function a_plus_equal_αb(a, α, b, n)
    # a, b same length, α float
    for i=1:n
        a[i] +=  α * b[i]
    end
    return a
end



function equilibrationRT(A, n_rows, n_cols)
    r = zeros(n_rows)
    t = zeros(n_cols)
    AT = sparse(A')
    for j=1:n_rows
        i = AT.colptr[j]
        k = AT.colptr[j+1] - 1
        r_j_inv = i <= k ? norm(AT.nzval[i:k], Inf) : 0.0
        if r_j_inv > 0.0
          r[j] = 1/r_j_inv
        end

    end

    for j=1:n_cols
        i = A.colptr[j]
        k = A.colptr[j+1] - 1
        t_j_inv = i <= k ? norm(A.nzval[i:k], Inf) : 0.0
        if t_j_inv > 0.0
          t[j] = 1/t_j_inv
        end
    end
    return Diagonal(r), Diagonal(t)
end


function mehrotraPCQuadBounds(QM; max_iter=100, ϵ_pdd=1e-8, ϵ_rb=1e-6, ϵ_rc=1e-6,
                              tol_Δx=1e-16, ϵ_μ=0., max_time=60., scaling=true,
                              display=true)

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
    n_rows, n_cols = size(A)
    @assert QM.meta.lcon == QM.meta.ucon # equality constraint (Ax=b)
    b = QM.meta.lcon
    Q = hess(QM, Oc)  # lower triangular

    if scaling
        R, T = equilibrationRT(A, n_rows, n_cols)
        inv_R_vec = zeros(n_rows)
        for i=1:n_rows
            if R[i,i] != 0.
                inv_R_vec[i] = 1/R[i,i]
            end
        end
        inv_T_vec = zeros(n_cols)
        for i=1:n_cols
            if T[i,i] != 0.
                inv_T_vec[i] = 1/T[i,i]
            end
        end
        inv_R = Diagonal(inv_R_vec)
        inv_T = Diagonal(inv_T_vec)
        A = R*A*T
        Q = Q*T^2
        c = T*c
        b = R*b
        lvar = inv_T*lvar
        uvar = inv_T*uvar
    end

    Arows, Acols, Avals = findnz(A)
    c0 = obj(QM, Oc)
    Qrows, Qcols, Qvals = findnz(Q)
    Q_sym = Symmetric(Q, :L)

    # init regularization values
    ρ, δ = 1e7*sqrt(eps()), 1e-1

    J_augmrows = vcat(Qcols, Acols, n_cols+1:n_cols+n_rows, 1:n_cols)
    J_augmcols = vcat(Qrows, Arows.+n_cols, n_cols+1:n_cols+n_rows, 1:n_cols)
    tmp_diag = zeros(n_cols)
    J_augmvals = vcat(-Qvals, Avals, δ*ones(n_rows), tmp_diag)
    J_augm = sparse(J_augmrows, J_augmcols, J_augmvals)

    x, λ, s_l, s_u, J_P = @views starting_points3(Q_sym, A, b, c, lvar, uvar, ilow, iupp,
                                     QM.meta.irng, J_augm, 1e-3, n_rows, n_cols)

#     x = init_x0_lsq(A, b, lvar, uvar)
#     @assert all(x .> lvar) && all(x .< uvar)
#     s_l, s_u = copy(Oc), copy(Oc)
#     s_l[ilow] .= 1.
#     s_u[iupp] .= 1.

    Qx = Q_sym * x
#     λ = sparse(A') \ (c+Qx) # least square initialisation, s_0 = stilde_0
    ATλ = A' * λ
    Ax = A * x
    rb = Ax - b
    rc = -Qx + ATλ + s_l - s_u - c

    x_minus_lvar = @views x[ilow] - lvar[ilow]
    uvar_minus_x = @views uvar[iupp] - x[iupp]
    μ = @views compute_μ(x_minus_lvar, uvar_minus_x,
                         s_l[ilow], s_u[iupp],
                         n_low, n_upp)

    k = 0
    Δ_aff = zeros(n_cols+n_rows+n_low+n_upp)
    Δ_cc = zeros(n_cols+n_rows+n_low+n_upp)
    Δ = zeros(n_cols+n_rows+n_low+n_upp)
    Δ_x_λ = zeros(n_cols+n_rows)

    # stopping criterion
    xTQx_2 = x' * Qx / 2.
    cTx = c' * x
    pri_obj = xTQx_2 + cTx + c0
    dual_obj = b' * λ - xTQx_2 + view(s_l,ilow)'*view(lvar,ilow) -
                    view(s_u,iupp)'*view(uvar,iupp) +c0
    pdd = abs(pri_obj - dual_obj ) / (1. + abs(pri_obj))
    max_rc, max_rb = norm(rc, Inf), norm(rb, Inf)
    optimal = pdd < ϵ_pdd && max_rb < ϵ_rb && max_rc < ϵ_rc

    n_Δx = 0.
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
        tmp_diag[ilow] .-= @views s_l[ilow] ./ x_minus_lvar
        tmp_diag[iupp] .-= @views s_u[iupp] ./ uvar_minus_x
        J_augmvals[end-n_cols+1:end] = tmp_diag

        J_augm = sparse(J_augmrows, J_augmcols, J_augmvals)

        #J_fact = ldlt(Symmetric(J_augm, :U))
        #J_fact = lu(Symmetric(J_augm), check=true)
        J_fact = ldl(Symmetric(J_augm, :U), J_P)

        Δ_aff = solve_augmented_system_aff!(J_fact, Δ_aff, Δ_x_λ, rc, rb, x_minus_lvar, uvar_minus_x,
                                            s_l, s_u, ilow, iupp,  n_cols, n_rows, n_low)

        α_aff_pri = @views compute_α_primal(x, Δ_aff[1:n_cols], lvar, uvar)
        α_aff_dual_l = @views compute_α_dual(s_l[ilow], Δ_aff[n_rows+n_cols+1: n_rows+n_cols+n_low])
        α_aff_dual_u = @views compute_α_dual(s_u[iupp], Δ_aff[n_rows+n_cols+n_low+1:end])

        # alpha_aff_dual_final is the min of the 2 alpha_aff_dual
        α_aff_dual_final = min(α_aff_dual_l, α_aff_dual_u)

        μ_aff = @views compute_μ(a_plus_αb(x_minus_lvar, α_aff_pri, Δ_aff[1:n_cols][ilow], n_low),
                                 a_plus_αb(uvar_minus_x, -α_aff_pri, Δ_aff[1:n_cols][iupp], n_upp),
                                 a_plus_αb(s_l[ilow], α_aff_dual_final, Δ_aff[n_rows+n_cols+1: n_rows+n_cols+n_low], n_low),
                                 a_plus_αb(s_u[iupp], α_aff_dual_final,  Δ_aff[n_rows+n_cols+n_low+1: end], n_upp),
                                 n_low, n_upp)

        σ = (μ_aff / μ)^3

        # corrector and centering step

        Δ_cc = solve_augmented_system_cc!(J_fact, Δ_cc, Δ_x_λ , Δ_aff, σ, μ,x_minus_lvar, uvar_minus_x,
                                          s_l, s_u, ilow, iupp, n_cols, n_rows, n_low)


        Δ = Δ_aff + Δ_cc # final direction

        α_pri = @views compute_α_primal(x, Δ[1:n_cols], lvar, uvar)
        α_dual_l = @views compute_α_dual(s_l[ilow], Δ[n_rows+n_cols+1: n_rows+n_cols+n_low])
        α_dual_u = @views compute_α_dual(s_u[iupp], Δ[n_rows+n_cols+n_low+1: end])

        α_dual_final = min(α_dual_l, α_dual_u)

        α_pri, α_dual_final = check_frontier!(0.05, x, x_minus_lvar, uvar_minus_x, s_l, s_u,
                                              Δ, α_pri, α_dual_final, lvar, uvar, ilow, iupp,
                                              n_low, n_upp, n_rows, n_cols)

        # new parameters
        x = @views a_plus_equal_αb(x, α_pri, Δ[1:n_cols], n_cols)
        λ = @views a_plus_equal_αb(λ, α_dual_final,
                                   Δ[n_cols+1: n_rows+n_cols], n_rows)
        s_l[ilow] = @views a_plus_equal_αb(s_l[ilow], α_dual_final,
                                           Δ[n_rows+n_cols+1: n_rows+n_cols+n_low], n_low)
        s_u[iupp] = @views a_plus_equal_αb(s_u[iupp], α_dual_final,
                                           Δ[n_rows+n_cols+n_low+1: end], n_upp)
        n_Δx = @views α_pri * norm(Δ[1:n_cols])
        x_minus_lvar = @views x[ilow] - lvar[ilow]
        uvar_minus_x = @views uvar[iupp] - x[iupp]

        μ = @views compute_μ(x_minus_lvar, uvar_minus_x,
                             s_l[ilow], s_u[iupp],
                             n_low, n_upp)

        Qx = mul!(Qx, Q_sym, x)
        xTQx_2 =  x' * Qx / 2.
        ATλ = mul!(ATλ, A', λ)
        Ax = mul!(Ax, A, x)
        cTx = c' * x
        pri_obj = xTQx_2 + cTx + c0
        dual_obj = b' * λ - xTQx_2 + view(s_l,ilow)'*view(lvar,ilow) -
                    view(s_u,iupp)'*view(uvar,iupp) +c0

        rb = Ax - b
        rc = -Qx + ATλ + s_l - s_u - c

        # update stopping criterion values:

        pdd = abs(pri_obj - dual_obj ) / (1. + abs(pri_obj))
        max_rc, max_rb = norm(rc, Inf), norm(rb, Inf)
        optimal = pdd < ϵ_pdd && max_rb < ϵ_rb && max_rc < ϵ_rc
        small_Δx, small_μ = n_Δx < tol_Δx, μ < ϵ_μ
        k += 1

        if δ >= 100*sqrt(eps())
            δ /= 10
            J_augmvals[end-n_cols-n_rows+1:end-n_cols] .= δ
        end
        if ρ >= 0.1*sqrt(eps())
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
        x = T*x
        Q_sym2 = Symmetric(Q*inv_T^2)
        Qx = mul!(Qx, Q_sym2, x)
        xTQx_2 =  x' * Qx / 2.
        A = rmul!(A, inv_T)
        A = lmul!(inv_R, A)
        b = mul!(b, inv_R, b)
        ATλ = mul!(ATλ, A', λ)
        Ax = mul!(Ax, A, x)
        cTx = c' *inv_T* x
        pri_obj = xTQx_2 + cTx + c0
        lvar = T * lvar
        uvar = T * uvar
        dual_obj = b' * λ - xTQx_2 + view(s_l,ilow)'*view(lvar,ilow) -
                    view(s_u,iupp)'*view(uvar,iupp) +c0
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
