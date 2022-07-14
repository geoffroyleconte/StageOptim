using QuadraticModels
Q = [6. 2. 1.
     2. 5. 2.
     1. 2. 4.]
c = [-8.; -3; -3]
A = [1. 0. 1.
     0. 2. 1.]
b = [0.; 3]
l = [0.; 0; 0]
u = [Inf; Inf; Inf]
QM = QuadraticModel(c, Q, A=A, lcon=b, ucon=b, lvar=l, uvar=u)

using RipQP
stats = ripqp(QM)

using LLSModels
m, n = 10, 15
A = rand(m, n)
b = rand(m)
lcon, ucon = zeros(m), fill(Inf, m)
C = ones(m, n)
lvar, uvar = fill(0.2, n), fill(1.0, n)
lls = LLSModel(A, b, lvar = lvar, uvar = uvar, C = C, lcon = lcon, ucon = ucon)
stats_lls = ripqp(lls)

path_pb = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\problemes_marosmeszaros"
path_to_QAFIRO = string(path_pb, "/QAFIRO.SIF")
using QPSReader
qps = readqps(path_to_QAFIRO)
qm = QuadraticModel(qps)
ripqp(qm)


stats = ripqp(qm, sp = K2KrylovParams(
  uplo = :U,
  kmethod = :minres_qlp,
  form_mat = true,
  equilibrate = true,
  preconditioner = LDL(T = Float32)),
)

ripqp(qm, solve_method = PC()) # default
ripqp(qm, solve_method = IPF())

path_pb_lp = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\problemes_netlib"
path_to_AGG = string(path_pb_lp, "/AGG.SIF")
qps2 = readqps(path_to_AGG)
qm2 = QuadraticModel(qps2)
using TimerOutputs
TimerOutputs.enable_debug_timings(RipQP)
reset_timer!(RipQP.to)
stats = ripqp(qm2, solve_method = PC(), sp = K2KrylovParams(
  uplo = :U,
  kmethod = :minres,
  form_mat = true,
  equilibrate = true,
  preconditioner = LDL(T = Float32)),
  display = false
)
TimerOutputs.complement!(RipQP.to)
show(RipQP.to, sortby = :firstexec)