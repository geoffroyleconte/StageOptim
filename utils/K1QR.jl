using RipQP, LinearAlgebra, LDLFactorizations, SparseArrays, QRMumps

include("./K1-utils.jl")

## K = [N̅⁻¹*fd.AT'; LH; M̅]'

struct K1QRParams{T<:Real} <: SolverParams
    method :: Symbol
    ρ_min  :: T
    δ_min  :: T
end

function K1QRParams(; method :: Symbol = :sne, ρ_min :: T = sqrt(eps()), δ_min :: T = sqrt(eps())) where {T<:Real} 
    if method != :sne && method != :ln 
        error("method should be :sne or :ln")  
    end 
    return K1QRParams{T}(method, ρ_min, δ_min) 
end

mutable struct PreallocatedData_K1QR{T<:Real} <: RipQP.PreallocatedData{T} 
    method      :: Symbol
    regu        :: RipQP.Regularization{T}
    K           :: SparseMatrixCSC{T, Int}
    Kqrm        :: qrm_spmat{T}
    K_fact      :: qrm_spfct{T}
    rhs         :: Vector{T}
    z           :: Vector{T}
    LH          :: SparseMatrixCSC{T, Int}
    D           :: Vector{T} 
    qp          :: Bool
end

function RipQP.PreallocatedData(sp :: K1QRParams, fd :: RipQP.QM_FloatData{T}, 
                                id :: RipQP.QM_IntData, 
                                iconf :: InputConfig{Tconf}) where {T<:Real, Tconf<:Real}

    D = one(T) .* ones(T, id.nvar)
    z = zeros(T, 2*id.nvar + id.ncon)
    rhs = zeros(T, id.nvar)
    regu = RipQP.Regularization(T(sqrt(eps())*1e5), T(sqrt(eps())*1e5), T(sp.ρ_min), T(sp.δ_min), :classic)

    qp = length(fd.Q.rowval) > 0
    if qp
        LDL = ldl(Symmetric(fd.Q, :U) + T(regu.ρ_min)*I)
        LH = LDL.L + I
        rmul!(LH, Diagonal(sqrt.(LDL.d)))
        permute!(LH, LDL.pinv, LDL.pinv)
        # TODO: create directly K without creating LDL.L+I and permute at the same time
    else 
        LH = spzeros(T, id.nvar, id.nvar)
    end

    K = create_K1_2(fd.AT, LH, D, regu.δ, id, qp)
    qrm_init()
    Kqrm = qrm_spmat_init(K)
    K_fact = qrm_spfct_init(Kqrm)
    sp.method == :sne && qrm_set(K_fact, "qrm_keeph", 0)
    qrm_analyse!(Kqrm, K_fact, transp='t')
    qrm_factorize!(Kqrm, K_fact, transp='t')

    return PreallocatedData_K1QR(sp.method, regu, K, Kqrm, K_fact, rhs, z, LH, D, qp)
end

function RipQP.update_pad!(pad :: PreallocatedData_K1QR{T}, dda :: RipQP.DescentDirectionAllocs{T}, 
                           pt :: RipQP.Point{T}, itd :: RipQP.IterData{T}, 
                           fd :: RipQP.Abstract_QM_FloatData{T}, id :: RipQP.QM_IntData, 
                           res :: RipQP.Residuals{T}, cnts :: RipQP.Counters, 
                           T0 :: RipQP.DataType) where {T<:Real}

    RipQP.update_regu!(pad.regu) 
    pad.D .= pad.regu.ρ
    pad.D[id.ilow] .+= pt.s_l ./ itd.x_m_lvar
    pad.D[id.iupp] .+= pt.s_u ./ itd.uvar_m_x
    update_K1_2!(pad.K.nzval, fd.AT.colptr, fd.AT.nzval, pad.D, pad.regu.δ, id.ncon, id.nvar)
    qrm_update!(pad.Kqrm, pad.K)
    qrm_factorize!(pad.Kqrm, pad.K_fact, transp='t')

end

function ln_qpsolve!(Δxy :: Vector{T}, K_fact :: qrm_spfct{T}, AT :: SparseMatrixCSC{T, I}, D :: Vector{T}, rhs :: Vector{T}, 
                     z :: Vector{T}, δ :: T, nvar :: Int, ncon :: Int, qp :: Bool) where {T <: AbstractFloat, I <: Integer}

    Δxy[nvar+1:end] ./= δ
    @views mul!(rhs, AT, Δxy[nvar+1:end]) 
    rhs .-= @views Δxy[1:nvar]

    qrm_solve!(K_fact, rhs, z, transp='t')
    qrm_apply!(K_fact, z)
    if qp
        Δxy[1:nvar] .= @views z[ncon+nvar+1:ncon+2*nvar] ./ sqrt.(D)
    else
        Δxy[1:nvar] .= @views z[ncon+1:ncon+nvar] ./ sqrt.(D)
    end
    @views mul!(Δxy[nvar+1:end], AT', Δxy[1:nvar], -one(T) / δ, one(T))
end

function sne_qpsolve!(Δxy :: Vector{T}, K_fact :: qrm_spfct{T}, AT :: SparseMatrixCSC{T, I}, D :: Vector{T}, rhs :: Vector{T}, 
                      z :: Vector{T}, δ :: T, nvar :: Int, ncon :: Int, qp :: Bool) where {T<:Real, I<:Integer}

    Δxy[nvar+1:end] ./= δ
    @views mul!(rhs, AT, Δxy[nvar+1:end]) 
    rhs .-= @views Δxy[1:nvar]
    qrm_solve!(K_fact, rhs, z, transp='t')
    qrm_solve!(K_fact, z, rhs)
    Δxy[1:nvar] .= rhs 
    @views mul!(Δxy[nvar+1:end], AT', Δxy[1:nvar], -one(T) / δ, one(T))
end


function RipQP.solver!(pad :: PreallocatedData_K1QR{T}, 
                       dda :: RipQP.DescentDirectionAllocsPC{T}, pt :: RipQP.Point{T}, 
                       itd :: RipQP.IterData{T}, fd :: RipQP.Abstract_QM_FloatData{T}, 
                       id :: RipQP.QM_IntData, res :: RipQP.Residuals{T}, 
                       cnts :: RipQP.Counters, T0 :: DataType, 
                       step :: Symbol) where {T<:Real}

    if step == :aff # affine predictor step
        if pad.method == :sne
            sne_qpsolve!(dda.Δxy_aff, pad.K_fact, fd.AT, pad.D, pad.rhs, pad.z, pad.regu.δ, id.nvar, id.ncon, pad.qp)
        else #:ln 
            ln_qpsolve!(dda.Δxy_aff, pad.K_fact, fd.AT, pad.D, pad.rhs, pad.z, pad.regu.δ, id.nvar, id.ncon, pad.qp)
        end
    else  
        if pad.method == :sne 
            sne_qpsolve!(itd.Δxy, pad.K_fact, fd.AT, pad.D, pad.rhs, pad.z, pad.regu.δ, id.nvar, id.ncon, pad.qp)
        else #:ln 
            ln_qpsolve!(itd.Δxy, pad.K_fact, fd.AT, pad.D, pad.rhs, pad.z, pad.regu.δ, id.nvar, id.ncon, pad.qp)
        end
    end

    return 0
end

using QuadraticModels, QPSReader
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
# qm = QuadraticModel(readqps(string(path_pb, "\\AGG.SIF"), mpsformat=:fixed))
# stats1 = ripqp(qm, iconf = InputConfig(sp = K1QRParams(method = :sne)))
# stats1bis = ripqp(qm, iconf = InputConfig(sp = K1QRParams(method = :ln)))
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
