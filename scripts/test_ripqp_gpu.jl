using QPSReader, QuadraticModels, SolverCore, CUDA, LinearAlgebra
CUDA.allowscalar(false)
path_pb = "/home/lecogeof/datasets/problemes_netlib"
# path_pb = "/home/lecogeof/datasets/problemes_netlib/problemes_marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"

using CUDA.CUSPARSE, SparseArrays
function QuadraticModelGPU(qps::QPSData; T = Float64)
  x0 = CUDA.zeros(T, qps.nvar)
  QuadraticModel(
    CuVector{T}(qps.c),
    CuSparseMatrixCSR{T}(sparse(qps.qrows, qps.qcols, qps.qvals, qps.nvar, qps.nvar)),
    A = CuSparseMatrixCSR{T}(sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)),
    lcon = CuVector{T}(qps.lcon),
    ucon = CuVector{T}(qps.ucon),
    lvar = CuVector{T}(qps.lvar),
    uvar = CuVector{T}(qps.uvar),
    c0 = T(qps.c0),
    x0 = x0,
  )
end
qm = QuadraticModelGPU(readqps(string(path_pb, "/AFIRO.SIF"), mpsformat=:fixed))
using RipQP
# include("/home/lecogeof/code/RipQP.jl/src/RipQP.jl")
stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                        sp = RipQP.K2StructuredParams(kmethod=:gpmr), 
                        solve_method = :IPF, presolve = false, scaling = false),
                     itol = RipQP.InputTol(max_iter=50, max_time=100.0,
                     ϵ_rc=1.0e-2, ϵ_rb=1.0e-2, ϵ_pdd=1.0e-2,))
# println(stats1)


# dense gpu
qm = QuadraticModel(readqps(string(path_pb, "/AFIRO.SIF"), mpsformat=:fixed))
using NLPModelsModifiers
qm = SlackModel(qm)
qm_gpu = QuadraticModel(
  CuVector{Float64}(qm.data.c),
  CUDA.CUSPARSE.CuSparseMatrixCSC(spzeros(Float64, qm.meta.nvar, qm.meta.nvar)),
  A = CUDA.CuMatrix{Float64}(qm.data.A),
  lcon = CuVector{Float64}(qm.meta.lcon),
  ucon = CuVector{Float64}(qm.meta.ucon),
  lvar = CuVector{Float64}(qm.meta.lvar),
  uvar = CuVector{Float64}(qm.meta.uvar),
  x0 = CUDA.zeros(Float64, qm.meta.nvar),
)

qm_gpu = QuadraticModel(
  CuVector{Float64}(qm.data.c),
  CUDA.CUSPARSE.CuSparseMatrixCSC(spzeros(Float64, qm.meta.nvar, qm.meta.nvar)),
  A = CUDA.CUSPARSE.CuSparseMatrixCSC(SparseMatrixCSC(qm.data.A)),
  lcon = CuVector{Float64}(qm.meta.lcon),
  ucon = CuVector{Float64}(qm.meta.ucon),
  lvar = CuVector{Float64}(qm.meta.lvar),
  uvar = CuVector{Float64}(qm.meta.uvar),
  x0 = CUDA.zeros(Float64, qm.meta.nvar),
)