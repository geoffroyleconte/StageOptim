mutable struct LDLLowPrecData{T <: Real, S, SI} <: PreconditionerDataK2{T, S}
  P::SI
end

function Identity(
  id::QM_IntData,
  fd::QM_FloatData{T},
  regu::Regularization{T},
  D::AbstractVector{T},
  K::LinearOperator{T},
) where {T <: Real}
  
  return IdentityData{T, typeof(fd.c), typeof(P)}(P)
end

function update_preconditioner!(
  pdat::IdentityData{T},
  pad::PreallocatedData{T},
  itd::IterData{T},
  pt::Point{T},
  id::QM_IntData,
  fd::QM_FloatData{T},
  cnts::Counters,
) where {T <: Real} end


struct K2LDLParams <: SolverParams
  regul::Symbol
end

function K2LDLParams(; regul::Symbol = :classic)
  regul == :classic ||
    regul == :dynamic ||
    regul == :none ||
    error("regul should be :classic or :dynamic or :none")
  return K2LDLParams(regul)
end

abstract type PreallocatedData_LDL{T <: Real, S} <: PreallocatedData{T, S} end

mutable struct PreallocatedData_K2LDL{T <: Real, S} <: PreallocatedData_LDL{T, S}
  D::S # temporary top-left diagonal
  regu::Regularization{T}
  diag_Q::SparseVector{T, Int} # Q diagonal
  K::SparseMatrixCSC{T, Int} # augmented matrix 
  K_fact::LDLFactorizations.LDLFactorization{T, Int, Int, Int} # factorized matrix
  fact_fail::Bool # true if factorization failed 
  diagind_K::Vector{Int} # diagonal indices of J
end

# outer constructor
function PreallocatedData(
  sp::K2LDLParams,
  fd::QM_FloatData{T},
  id::QM_IntData,
  iconf::InputConfig{Tconf},
) where {T <: Real, Tconf <: Real}

  # init Regularization values
  D = similar(fd.c, id.nvar)
  if iconf.mode == :mono
    regu = Regularization(
      T(sqrt(eps()) * 1e5),
      T(sqrt(eps()) * 1e5),
      1e-5 * sqrt(eps(T)),
      1e0 * sqrt(eps(T)),
      sp.regul,
    )
    D .= -T(1.0e0) / 2
  else
    regu = Regularization(
      T(sqrt(eps()) * 1e5),
      T(sqrt(eps()) * 1e5),
      T(sqrt(eps(T)) * 1e0),
      T(sqrt(eps(T)) * 1e0),
      sp.regul,
    )
    D .= -T(1.0e-2)
  end
  diag_Q = get_diag_Q(fd.Q.colptr, fd.Q.rowval, fd.Q.nzval, id.nvar)
  K = create_K2(id, D, fd.Q, fd.A, diag_Q, regu)

  diagind_K = get_diag_sparseCSC(K.colptr, id.ncon + id.nvar)
  K_fact = ldl_analyze(Symmetric(K, :U))
  if regu.regul == :dynamic
    Amax = @views norm(K.nzval[diagind_K], Inf)
    regu.ρ, regu.δ = -T(eps(T)^(3 / 4)), T(eps(T)^(0.45))
    K_fact.r1, K_fact.r2 = regu.ρ, regu.δ
    K_fact.tol = Amax * T(eps(T))
    K_fact.n_d = id.nvar
  elseif regu.regul == :none
    regu.ρ, regu.δ = zero(T), zero(T)
  end
  K_fact = ldl_factorize!(Symmetric(K, :U), K_fact)
  K_fact.__factorized = true

  return PreallocatedData_K2LDL(
    D,
    regu,
    diag_Q, #diag_Q
    K, #K
    K_fact, #K_fact
    false,
    diagind_K, #diagind_K
  )
end

function convertpad(
  ::Type{<:PreallocatedData{T}},
  pad::PreallocatedData_K2LDL{T_old},
  T0::DataType,
) where {T <: Real, T_old <: Real}
  pad = PreallocatedData_K2LDL(
    convert(Array{T}, pad.D),
    convert(Regularization{T}, pad.regu),
    convert(SparseVector{T, Int}, pad.diag_Q),
    convert(SparseMatrixCSC{T, Int}, pad.K),
    convertldl(T, pad.K_fact),
    pad.fact_fail,
    pad.diagind_K,
  )

  if pad.regu.regul == :classic
    if T == Float64 && T0 == Float64
      pad.regu.ρ_min, pad.regu.δ_min = T(sqrt(eps()) * 1e-5), T(sqrt(eps()) * 1e0)
    else
      pad.regu.ρ_min, pad.regu.δ_min = T(sqrt(eps(T)) * 1e1), T(sqrt(eps(T)) * 1e1)
    end
  elseif pad.regu.regul == :dynamic
    pad.regu.ρ, pad.regu.δ = -T(eps(T)^(3 / 4)), T(eps(T)^(0.45))
    pad.K_fact.r1, pad.K_fact.r2 = pad.regu.ρ, pad.regu.δ
  end

  return pad
end

# function used to solve problems
# solver LDLFactorization
function solver!(
  pad::PreallocatedData_K2LDL{T},
  dda::DescentDirectionAllocs{T},
  pt::Point{T},
  itd::IterData{T},
  fd::Abstract_QM_FloatData{T},
  id::QM_IntData,
  res::AbstractResiduals{T},
  cnts::Counters,
  T0::DataType,
  step::Symbol,
) where {T <: Real}
  # erase dda.Δxy_aff only for affine predictor step with PC method
  step == :aff ? ldiv!(pad.K_fact, dda.Δxy_aff) : ldiv!(pad.K_fact, itd.Δxy)
  return 0
end

function update_pad!(
  pad::PreallocatedData_K2LDL{T},
  dda::DescentDirectionAllocs{T},
  pt::Point{T},
  itd::IterData{T},
  fd::Abstract_QM_FloatData{T},
  id::QM_IntData,
  res::AbstractResiduals{T},
  cnts::Counters,
  T0::DataType,
) where {T <: Real}
  if pad.regu.regul == :classic && cnts.k != 0 # update ρ and δ values, check K diag magnitude 
    out = update_regu_diagK2!(
      pad.regu,
      pad.K.nzval,
      pad.diagind_K,
      id.nvar,
      itd.pdd,
      itd.l_pdd,
      itd.mean_pdd,
      cnts,
      T,
      T0,
    )
    out == 1 && return out
  end

  out = factorize_K2!(
    pad.K,
    pad.K_fact,
    pad.D,
    pad.diag_Q,
    pad.diagind_K,
    pad.regu,
    pt.s_l,
    pt.s_u,
    itd.x_m_lvar,
    itd.uvar_m_x,
    id.ilow,
    id.iupp,
    id.ncon,
    id.nvar,
    cnts,
    itd.qp,
    T,
    T0,
  ) # update D and factorize K

  if out == 1
    pad.fact_fail = true
    return out
  end

  return 0
end

