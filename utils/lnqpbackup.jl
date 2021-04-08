
using RipQP, LinearAlgebra, LDLFactorizations, SparseArrays, qr_mumps

function ln_qp(K_fact :: qrm_spfct{T}, z :: Vector{T}, A , D :: Diagonal{T}, N :: Diagonal{T}, ξp :: Vector{T}, ξd :: Vector{T}) where {T <: AbstractFloat, I <: Integer}
    m, n = size(A)
    Aᵀ = A'
    # M⁻¹ = inv(M)
    N⁻¹ = inv(N)
    M̅ = Diagonal(sqrt.(D.diag))
    N̅ = Diagonal(sqrt.(N.diag))
    M̅⁻¹ = inv(M̅)
    N̅⁻¹ = inv(N̅)
    # z = zeros(T, 2*n + m)
    # F = qr([N̅⁻¹*A; LHᵀ; M̅]) 
    # L = LowerTriangular(F.R')
    rhs = Aᵀ * (N⁻¹ * ξp) - ξd
    # permute!(rhs, F.pcol)
    # z[1:n] = L \ rhs
    z .= 0.
    qrm_solve!(K_fact, rhs, z, transp='t')
    # Δs = F.Q * z
    Δs = qrm_apply(K_fact, z)
    # invpermute!(Δs, F.prow)
    Δx = M̅⁻¹ * Δs[m+n+1:m+2*n]
    Δy = N⁻¹ * (ξp - A * Δx)
    return Δx, Δy
end

function ln_qpsolve!(Δxy :: Vector{T}, K_fact :: qrm_spfct{T}, AT :: SparseMatrixCSC{T, I}, D :: Vector{T}, rhs :: Vector{T}, 
                     z :: Vector{T}, δ :: T, nvar :: Int, ncon :: Int) where {T <: AbstractFloat, I <: Integer}
    # M⁻¹ = inv(M)
    # N⁻¹ = inv(N)
    # M̅ = Diagonal(sqrt.(D.diag))
    # N̅ = Diagonal(sqrt.(N.diag))
    # M̅⁻¹ = inv(M̅)
    # N̅⁻¹ = inv(N̅)
    # z = zeros(T, 2*n + m)
    # F = qr([N̅⁻¹*A; LHᵀ; M̅]) 
    # L = LowerTriangular(F.R')
    Δxy[nvar+1:end] ./= δ
    # rhs = Aᵀ * (N⁻¹ * ξp) - ξd
    @views mul!(rhs, AT, Δxy[nvar+1:end]) 
    rhs .-= @views Δxy[1:nvar]
    # permute!(rhs, F.pcol)
    # z[1:n] = L \ rhs
    # z .= 0.
    qrm_solve!(K_fact, rhs, z, transp='t')
    # Δs = F.Q * z
    qrm_apply!(K_fact, z)
    # invpermute!(Δs, F.prow)
    # Δx = M̅⁻¹ * Δs[m+n+1:m+2*n]
    # Δy = N⁻¹ * (ξp - A * Δx)
    Δxy[1:nvar] .= @views z[ncon+nvar+1:ncon+2*nvar] ./ sqrt.(D)
    @views mul!(Δxy[nvar+1:end], AT', Δxy[1:nvar], -one(T) / δ, one(T))
end

struct K1LNQPParams{T<:Real} <: SolverParams
    variant :: Bool
    ρ :: T # dual regularization 
    δ :: T # primal regularization
end

mutable struct PreallocatedData_K1LNQP{T<:Real} <: RipQP.PreallocatedData{T} 
    sp          :: K1LNQPParams{T}
    K           :: SparseMatrixCSC{T, Int}
    Kqrm        :: qrm_spmat{T}
    K_fact      :: qrm_spfct{T}
    rhs         :: Vector{T}
    z           :: Vector{T}
    LH          :: SparseMatrixCSC{T, Int}
    D           :: Vector{T} 
    N           :: Vector{T}
    ξp          :: Vector{T}
    ξd          :: Vector{T}
end

function create_K1_2(AT, LH, dH, pinvH, D, δ, id)
    T = eltype(D)
    nnzK = length(AT.rowval) + length(LH.rowval) + id.nvar # A, LH, D
    K_colptr = Vector{Int}(undef, id.ncon+2*id.nvar+1) 
    K_rowval = Vector{Int}(undef, nnzK) 
    K_nzval = Vector{T}(undef, nnzK) 

    fill_K1_2!(K_colptr, K_rowval, K_nzval, AT.colptr, AT.rowval, AT.nzval, LH.colptr, LH.rowval, LH.nzval, 
               dH, pinvH, D, δ, id.nvar, id.ncon)

    return SparseMatrixCSC(id.nvar, id.ncon+2*id.nvar, K_colptr, K_rowval, K_nzval)
end

function fill_K1_2!(K_colptr, K_rowval, K_nzval, AT_colptr, AT_rowval, AT_nzval, LH_colptr, LH_rowval, LH_nzval, dH, pinvH, 
                    D, δ, nvar, ncon)

    # fill first bloc
    K_colptr[1:ncon+1] = AT_colptr
    sqrtδ = sqrt(δ)
    for j=1:ncon
        for k=AT_colptr[j]:(AT_colptr[j+1]-1)
            K_rowval[k] = AT_rowval[k]
            K_nzval[k] = AT_nzval[k] / sqrtδ
        end
    end

    # 2nd bloc
    K_deb = K_colptr[ncon+1]
    for i=1:nvar
        K_colptr[ncon+1+i] = K_deb + LH_colptr[i+1] - 1
    end
    for j=1:nvar 
        for k=LH_colptr[j]: (LH_colptr[j+1]-1)
            idx_K = K_deb + k - 1
            K_rowval[idx_K] = LH_rowval[k]
            K_nzval[idx_K] = LH_nzval[k]
        end
    end

    # third bloc
    K_deb = K_colptr[ncon+nvar+1]
    for i=1:nvar
        idx_K = K_deb + i 
        K_colptr[ncon+nvar+i+1] = idx_K
        K_rowval[idx_K-1] = i 
        K_nzval[idx_K-1] = sqrt(D[i])
    end
end

function RipQP.PreallocatedData(sp :: SolverParams, fd :: RipQP.QM_FloatData{T}, 
                                id :: RipQP.QM_IntData, 
                                iconf :: InputConfig{Tconf}) where {T<:Real, Tconf<:Real}

    D = one(T) .* ones(T, id.nvar)
    N = sp.δ .* ones(T, id.ncon)
    ξp, ξd = copy(fd.b), zeros(T, id.nvar)
    z = zeros(T, 2*id.nvar + id.ncon)
    rhs = zeros(T, id.nvar)

    LDL = ldl(Symmetric(fd.Q, :U) + T(sp.ρ)*I)
    LH = LDL.L + I
    rmul!(LH, Diagonal(sqrt.(LDL.d)))
    permute!(LH, LDL.pinv, LDL.pinv)

    # form K1
    N⁻¹ = Diagonal(one(T) ./ N)
    M̅ = Diagonal(sqrt.(D))
    N̅ = Diagonal(sqrt.(N))
    M̅⁻¹ = inv(M̅)
    N̅⁻¹ = inv(N̅)
    # K = [N̅⁻¹*fd.AT'; LH; M̅]
    # K = hcat(fd.AT*N̅⁻¹, LH, M̅)
    K = create_K1_2(fd.AT, LH, LDL.d, LDL.pinv, D, sp.δ, id)

    qrm_init()
    Kqrm = qrm_spmat_init(K)
    K_fact = qrm_spfct_init(Kqrm)
    qrm_analyse!(Kqrm, K_fact, transp='t')
    qrm_factorize!(Kqrm, K_fact, transp='t')

    return PreallocatedData_K1LNQP(sp, K, Kqrm, K_fact, rhs, z, LH, D, N, ξp, ξd)
end

function update_K1_2!(K_nzval, AT_colptr, AT_nzval, D, δ, ncon, nvar)

    sqrtδ = sqrt(δ)
    for j=1:ncon
        for k=AT_colptr[j]:(AT_colptr[j+1]-1)
            K_nzval[k] = AT_nzval[k] / sqrtδ
        end
    end

    K_nzval[end-nvar+1:end] .= sqrt.(D)
end

function RipQP.update_pad!(pad :: PreallocatedData_K1LNQP{T}, dda :: RipQP.DescentDirectionAllocs{T}, 
                           pt :: RipQP.Point{T}, itd :: RipQP.IterData{T}, 
                           fd :: RipQP.Abstract_QM_FloatData{T}, id :: RipQP.QM_IntData, 
                           res :: RipQP.Residuals{T}, cnts :: RipQP.Counters, 
                           T0 :: RipQP.DataType) where {T<:Real}

    pad.D .= pad.sp.ρ
    pad.D[id.ilow] .+= pt.s_l ./ itd.x_m_lvar
    pad.D[id.iupp] .+= pt.s_u ./ itd.uvar_m_x
    N⁻¹ = Diagonal(one(T) ./ pad.N)
    M̅ = Diagonal(sqrt.(pad.D))
    N̅ = Diagonal(sqrt.(pad.N))
    M̅⁻¹ = inv(M̅)
    N̅⁻¹ = inv(N̅)
    # pad.K = [N̅⁻¹*fd.AT'; pad.LH; M̅]
    # pad.K = hcat(fd.AT*N̅⁻¹, pad.LH, M̅)
    update_K1_2!(pad.K.nzval, fd.AT.colptr, fd.AT.nzval, pad.D, pad.sp.δ, id.ncon, id.nvar)
    qrm_update!(pad.Kqrm, pad.K)
    qrm_factorize!(pad.Kqrm, pad.K_fact, transp='t')

end

function RipQP.solver!(pad :: PreallocatedData_K1LNQP{T}, 
                       dda :: RipQP.DescentDirectionAllocsPC{T}, pt :: RipQP.Point{T}, 
                       itd :: RipQP.IterData{T}, fd :: RipQP.Abstract_QM_FloatData{T}, 
                       id :: RipQP.QM_IntData, res :: RipQP.Residuals{T}, 
                       cnts :: RipQP.Counters, T0 :: DataType, 
                       step :: Symbol) where {T<:Real}

    if step == :aff # affine predictor step
        # pad.ξp .= @views dda.Δxy_aff[id.nvar+1:end]
        # pad.ξd .= @views dda.Δxy_aff[1:id.nvar]
 
        # # Δx_aff, Δy_aff = ln_qp(fd.AT', pad.LH, Diagonal(pad.D), Diagonal(pad.N), pad.ξp, pad.ξd)
        # Δx_aff, Δy_aff = ln_qp(pad.K_fact, pad.z, fd.AT', Diagonal(pad.D), Diagonal(pad.N), pad.ξp, pad.ξd)
        # dda.Δxy_aff[1:id.nvar] .= Δx_aff
        # dda.Δxy_aff[id.nvar+1:end] .= Δy_aff
        ln_qpsolve!(dda.Δxy_aff, pad.K_fact, fd.AT, pad.D, pad.rhs, pad.z, pad.sp.δ, id.nvar, id.ncon)
    else  
        # pad.ξp .= @views itd.Δxy[id.nvar+1:end]
        # pad.ξd .= @views itd.Δxy[1:id.nvar]
        # # Δx, Δy = ln_qp(fd.AT', pad.LH, Diagonal(pad.D), Diagonal(pad.N), pad.ξp, pad.ξd)
        # Δx, Δy = ln_qp(pad.K_fact, pad.z, fd.AT', Diagonal(pad.D), Diagonal(pad.N), pad.ξp, pad.ξd)
        # itd.Δxy[1:id.nvar] .= Δx
        # itd.Δxy[id.nvar+1:end] .= Δy
        ln_qpsolve!(itd.Δxy, pad.K_fact, fd.AT, pad.D, pad.rhs, pad.z, pad.sp.δ, id.nvar, id.ncon)
    end

    return 0
end

using QuadraticModels, QPSReader
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
qm = QuadraticModel(readqps(string(path_pb, "\\QSEBA.SIF"), mpsformat=:fixed))
stats1 = ripqp(qm, iconf = InputConfig(sp = K1LNQPParams(true, 1.0e-8, 1.0e-8)))
# stats2 = ripqp(qm)

# Q = [6. 2. 1.
#      2. 5. 2.
#      1. 2. 4.]
# c = [-8.; -3; -3]
# A = [1. 0. 1.
#      0. 2. 1.]
# b = [0.; 3]
# l = [0.;0;0]
# u = [Inf; Inf; Inf]
# QM = QuadraticModel(c, Q, A=A, lcon=b, ucon=b, lvar=l, uvar=u, c0=0., name="QM1");
# stats2 = ripqp(QM, iconf = InputConfig(sp = K1LNQPParams(true, 1.0e-8, 1.0e-8)))
