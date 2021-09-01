# a Basic LDLᵀ solver
using RipQP, LinearAlgebra, LDLFactorizations, SparseArrays

struct K2basicLDLParams{T<:Real} <: SolverParams
    ρ :: T
    δ :: T
end

mutable struct PreallocatedData_K2basic{T<:Real} <: RipQP.PreallocatedData{T} 
    D                :: Vector{T}                                        # temporary top-left diagonal
    ρ                :: T
    δ                :: T
    K                :: SparseMatrixCSC{T,Int} # augmented matrix 
    K_fact           :: LDLFactorizations.LDLFactorization{T,Int,Int,Int} # factorized matrix
end

# outer constructor
function RipQP.PreallocatedData(sp :: SolverParams, fd :: RipQP.QM_FloatData{T}, id :: RipQP.QM_IntData, 
                                iconf :: InputConfig{Tconf}) where {T<:Real, Tconf<:Real}

    ρ, δ = T(sp.ρ), T(sp.δ)
    K = spzeros(T, id.ncon+id.nvar, id.ncon + id.nvar)
    K[1:id.nvar, 1:id.nvar] = .-fd.Q .- ρ .* Diagonal(ones(T, id.nvar))     
    K[1:id.nvar, id.nvar+1:end] = fd.AT      
    K[diagind(K)[id.nvar+1:end]] .= δ

    K_fact = ldl_analyze(Symmetric(K, :U))
    K_fact = ldl_factorize!(Symmetric(K, :U), K_fact)
    K_fact.__factorized = true

    return PreallocatedData_K2basic(zeros(T, id.nvar),
                                    ρ,
                                    δ,
                                    K, #K
                                    K_fact #K_fact
                                    )
end

function RipQP.convertpad(::Type{<:RipQP.PreallocatedData{T}}, pad :: PreallocatedData_K2basic{T_old}, 
                          T0 :: DataType) where {T<:Real, T_old<:Real} 

    pad = PreallocatedData_K2basic(convert(Array{T}, pad.D),
                                   T(pad.ρ),
                                   T(pad.δ),
                                   convert(SparseMatrixCSC{T,Int}, pad.K),
                                   RipQP.convertldl(T, pad.K_fact)
                                   )

    pad.ρ /= 100
    pad.δ /= 100
    return pad
end

function RipQP.update_pad!(pad :: PreallocatedData_K2basic{T}, dda :: RipQP.DescentDirectionAllocs{T}, 
                           pt :: RipQP.Point{T}, itd :: RipQP.IterData{T}, fd :: RipQP.Abstract_QM_FloatData{T}, 
                           id :: RipQP.QM_IntData, res :: RipQP.Residuals{T}, cnts :: RipQP.Counters, 
                           T0 :: RipQP.DataType) where {T<:Real}

    # update the diagonal of K2
    pad.D .= -pad.ρ
    pad.D[id.ilow] .-= pt.s_l ./ itd.x_m_lvar
    pad.D[id.iupp] .-= pt.s_u ./ itd.uvar_m_x
    pad.D .-= fd.Q[diagind(fd.Q)]
    pad.K[diagind(pad.K)[1:id.nvar]] = pad.D 
    pad.K[diagind(pad.K)[id.nvar+1:end]] .= pad.δ

    # factorize K2
    ldl_factorize!(Symmetric(pad.K, :U), pad.K_fact)

end

# function used to solve problems
# solver LDLFactorization
function RipQP.solver!(pad :: PreallocatedData_K2basic{T}, 
                       dda :: RipQP.DescentDirectionAllocsPC{T}, pt :: RipQP.Point{T}, 
                       itd :: RipQP.IterData{T}, fd :: RipQP.Abstract_QM_FloatData{T}, 
                       id :: RipQP.QM_IntData, res :: RipQP.Residuals{T}, 
                       cnts :: RipQP.Counters, T0 :: DataType, 
                       step :: Symbol) where {T<:Real}
    
    if step == :aff # affine predictor step
        # solve the system and overwrite dda.Δxy_aff
        ldiv!(pad.K_fact, dda.Δxy_aff) 
    else # for all other steps including the initial point
        # solve the system and overwrite itd.Δxy
        ldiv!(pad.K_fact, itd.Δxy)
    end

    return 0
end

using QuadraticModels, QPSReader
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
qm = QuadraticModel(readqps(string(path_pb, "\\QSEBA.SIF"), mpsformat=:fixed))
stats1 = ripqp(qm, iconf = RipQP.InputConfig(refinement = :none, kc=0,mode=:multi, scaling=true,
                    sp = K2basicLDLParams(1.0e-6, 1.0e-6), solve_method=:PC),
                     itol = RipQP.InputTol(max_iter=100, ϵ_rb32 = 1e-6) )#