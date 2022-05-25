using QuadraticModels, QPSReader, SolverTools
using Quadmath, SparseArrays
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\quadLP\\data\\MPS"
path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"

function createQuadraticModel128(qpdata; name="qp_pb")
    return QuadraticModel(convert(Array{Float128}, qps1.c), qpdata.qrows, qpdata.qcols,
            convert(Array{Float128}, qps1.qvals),
            Arows=qpdata.arows, Acols=qpdata.acols,
            Avals=convert(Array{Float128}, qps1.avals),
            lcon=convert(Array{Float128}, qps1.lcon),
            ucon=convert(Array{Float128}, qps1.ucon),
            lvar=convert(Array{Float128}, qps1.lvar),
            uvar=convert(Array{Float128}, qps1.uvar),
            c0=Float128(qpdata.c0), x0 = zeros(Float128, length(qps1.c)), name=name)
end

qps1 = readqps(string(path_pb, "/TMA_ME_presolved.mps"))
# qps1 = readqps(string(path_pb, "/GlcAerWT_presolved.mps"))
# qps1 = readqps(string(path_pb, "/GlcAlift.mps"))

qm1 = createQuadraticModel128(qps1)

using RipQP
stats1 = RipQP.ripqp(qm1, 
  mode = :multiref,
  Timulti = Float64,
  sp = K2LDLParams(regul = :hybrid, ρ_min=1.0e-7, δ_min = 1.0e-6),
  sp3 = RipQP.K2KrylovParams{Float128}(
    uplo = :U,
    kmethod=:gmres,
    rhs_scale=true, #δ0 = 1.0e-20,ρ0=1.0e-20,
    form_mat = true,
    equilibrate = false,
    mem = 100,
    preconditioner = LDLLowPrec(T = Float64, warm_start = true, pos = :L),
    ρ_min=Float128(1.0e-16),
    δ_min = Float128(1.0e-19),
    itmax = 0,
    atol_min = Float128(1.0e-16),
    rtol_min = Float128(1.0e-16),
  ),
  solve_method=RipQP.IPF(γ = 0.1, r = 0.95),
  scaling = true,
  history=false,
  perturb = true,
  ps=false,
  itol = RipQP.InputTol(
    Float128,
    ϵ_rb = Float128(1e-40),
    max_iter = 1500,
    max_time = 70000.0,
    max_iter64 = 200,
  ),
  display = true,
)

println(stats1)

