using QPSReader, QuadraticModels, SolverCore, CUDA, LinearAlgebra, DoubleFloats
CUDA.allowscalar(false)
path_pb = "/home/lecogeof/datasets/problemes_netlib"
# path_pb = "/home/lecogeof/datasets/problemes_netlib/problemes_marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"

using CUDA.CUSPARSE, SparseArrays
function QuadraticModelGPU(qps::QPSData; T = Float64)
  x0 = CUDA.zeros(T, qps.nvar)
  QuadraticModel(
    CuVector{T}(qps.c),
    CuSparseMatrixCSC{T}(sparse(qps.qrows, qps.qcols, qps.qvals, qps.nvar, qps.nvar)),
    A = CuSparseMatrixCSC{T}(sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)),
    lcon = CuVector{T}(qps.lcon),
    ucon = CuVector{T}(qps.ucon),
    lvar = CuVector{T}(qps.lvar),
    uvar = CuVector{T}(qps.uvar),
    c0 = T(qps.c0),
    x0 = x0,
  )
end
function QuadraticModelGPUCOO(qps::QPSData; T = Float64)
  x0 = CUDA.zeros(T, qps.nvar)
  Ti = (T == Float32) ? Int32 : Int64
  QuadraticModel(
    CuVector{T}(qps.c),
    CuSparseMatrixCOO{T, Ti}(CuVector(qps.qrows), CuVector(qps.qcols), CuVector(qps.qvals), (qps.nvar, qps.nvar), length(qps.qrows)),
    A = CuSparseMatrixCOO{T, Ti}(CuVector(qps.arows), CuVector(qps.acols), CuVector(qps.avals), (qps.ncon, qps.nvar), length(qps.arows)),
    lcon = CuVector{T}(qps.lcon),
    ucon = CuVector{T}(qps.ucon),
    lvar = CuVector{T}(qps.lvar),
    uvar = CuVector{T}(qps.uvar),
    c0 = T(qps.c0),
    x0 = x0,
  )
end
# qm = QuadraticModelGPU(readqps(string(path_pb, "/AFIRO.SIF"), mpsformat=:fixed))
T = Double64
qm = QuadraticModelGPUCOO(readqps(string(path_pb, "/AGG.SIF"), mpsformat=:fixed), T = T)
using RipQP
# # include("/home/lecogeof/code/RipQP.jl/src/RipQP.jl")
stats1 = RipQP.ripqp(qm, 
                     sp = RipQP.K2KrylovGPUParams(kmethod=:minres, uplo = :U, preconditioner = RipQP.LDLGPU(T = Float32), verbose = 0), 
                     solve_method = IPF(), ps = false, scaling = false,
                     itol = RipQP.InputTol(T, max_iter=50, max_time=100.0, ϵ_rc=T(1.0e-4), ϵ_rb=T(1.0e-4), ϵ_pdd=T(1.0e-4)))
# println(stats1)
# stats1 = RipQP.ripqp(qm, 
#                      sp = RipQP.K2KrylovGPUParams(kmethod=:minres, uplo = :U, equilibrate = false), 
#                      solve_method = IPF(), ps = false, scaling = false,
#                      itol = RipQP.InputTol(max_iter=50, max_time=100.0, ϵ_rc=1.0e-2, ϵ_rb=1.0e-2, ϵ_pdd=1.0e-2))

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



using QPSReader, QuadraticModels, SolverCore, CUDA, LinearAlgebra
CUDA.allowscalar(false)
path_pb = "/home/lecogeof/datasets/problemes_netlib"
# path_pb = "/home/lecogeof/datasets/problemes_netlib/problemes_marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"

using CUDA.CUSPARSE, SparseArrays
Q = [6. 2. 1.
     2. 5. 2.
     1. 2. 4.]
c = [-8.; -3; -3]
A = [1. 0. 1.
     0. 2. 1.]
b = [0.; 3]
l = [0.;0;0]
u = [Inf; Inf; Inf]
Qsp = sparse(tril(Q))
Qrows, Qcols, Qvals = findnz(Qsp)
Asp = sparse(A)
Arows, Acols, Avals = findnz(Asp)
x0 = CUDA.zeros(Float64, 3)
QM = QuadraticModel(
  CuVector{Float64}(c),
  CuSparseMatrixCOO{Float64, Int64}(CuVector(Qrows), CuVector(Qcols), CuVector(Qvals), size(Q), nnz(Qsp)),
  A = CuSparseMatrixCOO{Float64, Int64}(CuVector(Arows), CuVector(Acols), CuVector(Avals), size(A), nnz(A)),
  lcon=CuVector{Float64}(b),
  ucon=CuVector{Float64}(b),
  lvar=CuVector{Float64}(l),
  uvar=CuVector{Float64}(u),
  x0 = x0,
);
using RipQP
stats1 = RipQP.ripqp(QM, 
                     sp = RipQP.K2KrylovGPUParams(kmethod=:gmres, uplo = :U, equilibrate = false), 
                     solve_method = IPF(), ps = false, scaling = false,
                     itol = RipQP.InputTol(max_iter=50, max_time=100.0, ϵ_rc=1.0e-2, ϵ_rb=1.0e-2, ϵ_pdd=1.0e-2))
# println(stats1)

A= [-15.8756  -2.0       -1.0     1.0         0.0;
-2.0     -9.82248   -2.0     0.0         2.0;
-1.0     -2.0      -11.8426  1.0         1.0;
 1.0      0.0        1.0     0.00149012  0.0;
 0.0      2.0        1.0     0.0         0.00149012]
b = [0.0; 0; 0; 0; 1.0]
Asp = Symmetric(sparse(triu(A)), :U)
(x, stats) = minres(Asp, b, verbose = 1)
println(Asp * x - b)
Agpu = Symmetric(CUDA.CUSPARSE.CuSparseMatrixCSC(Asp.data), :U);
bgpu = CuVector{Float64}(b)
(xgpu, stats) = minres(Agpu, bgpu, verbose = 1)
println(Agpu * xgpu - bgpu)