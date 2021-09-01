using QPSReader, QuadraticModels, SolverCore, CUDA
CUDA.allowscalar(false)
include("/home/lecogeof/code/RipQP/src/RipQP.jl")
path_pb = "/home/lecogeof/datasets/problemes_netlib"
# path_pb = "/home/lecogeof/datasets/problemes_netlib/problemes_marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"

function QuadraticModelGPU(qps::QPSData, x0 = CUDA.zeros(Float64, qps.nvar))
  QuadraticModel(
    CuVector(qps.c),
    qps.qrows,
    qps.qcols,
    CuVector(qps.qvals),
    Arows = qps.arows,
    Acols = qps.acols,
    Avals = CuVector(qps.avals),
    lcon = CuVector(qps.lcon),
    ucon = CuVector(qps.ucon),
    lvar = CuVector(qps.lvar),
    uvar = CuVector(qps.uvar),
    c0 = qps.c0,
    x0 = x0,
  )
end
qm = QuadraticModelGPU(readqps(string(path_pb, "/AFIRO.SIF"), mpsformat=:fixed))
stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                        sp = RipQP.K2minresParams(preconditioner = :Jacobi, ratol=1.0e-6, rrtol=1.0e-6), scaling = false),
                     itol = RipQP.InputTol(max_iter=50, max_time=20.0,
                     ϵ_rc=1.0e-2, ϵ_rb=1.0e-2, ϵ_pdd=1.0e-2,))
# println(stats1)