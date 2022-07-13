using QuadraticModels, QPSReader, SolverTools
using Quadmath, SparseArrays, DoubleFloats
path_pb = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\quadLP\data\MPS"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"

function createQuadraticModel_T(qpdata; T = Float128, name="qp_pb")
    return QuadraticModel(convert(Array{T}, qps1.c), qpdata.qrows, qpdata.qcols,
            convert(Array{T}, qps1.qvals),
            Arows=qpdata.arows, Acols=qpdata.acols,
            Avals=convert(Array{T}, qps1.avals),
            lcon=convert(Array{T}, qps1.lcon),
            ucon=convert(Array{T}, qps1.ucon),
            lvar=convert(Array{T}, qps1.lvar),
            uvar=convert(Array{T}, qps1.uvar),
            c0=T(qpdata.c0), x0 = zeros(T, length(qps1.c)), name=name)
end

# qps1 = readqps(string(path_pb, "CYCLE.SIF"))
qps1 = readqps(string(path_pb, "/TMA_ME_presolved.mps"))
# qps1 = readqps(string(path_pb, "/GlcAerWT_presolved.mps"))
# qps1 = readqps(string(path_pb, "/GlcAlift.mps"))

using RipQP
# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
T = Float128
qm1 = createQuadraticModel_T(qps1, T = T) # create QuadraticModel in Float128
Tlow = Float64
stats1 = ripqp(qm1, 
  mode = :multi,
  Timulti = Tlow,
  sp = K2LDLParams(regul = :hybrid, ρ_min=1.0e-8, δ_min = 1.0e-8), # solve in Float64
  sp3 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:minres_qlp,
    form_mat = true,
    equilibrate = true,
    itmax = 10,
    preconditioner = LDL(T = Tlow, warm_start = true, pos = :L),
    ρ_min=T(1.0e-16),
    δ_min = T(1.0e-19),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(γ = 0.1, r = 0.95),
  ps=false,
  itol = InputTol(
    T,
    ϵ_rb = T(1e-40), # very small to see what residuals can be reached
    max_iter = 300,
    max_time = 70000.0,
    max_iter64 = 200,
  ),
  display = true,
)

println(stats1)

