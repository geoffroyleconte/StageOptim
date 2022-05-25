using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, SparseArrays, TimerOutputs, PrettyTables
# using RipQP
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib_ps"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
TimerOutputs.enable_debug_timings(RipQP)
max_time = 50.0
max_iter = 100
kmethod = :dqgmres
ftype = Float32
ρ_min, δ_min = 1.0e-7, 1.0e-8
atol0, rtol0 = 1.0e-2, 1.0e-2
atol_min, rtol_min = 1.0e-4, 1.0e-4

# compile
qm = QuadraticModel(readqps(string(path_pb, "\\", "AFIRO", ".SIF"), mpsformat=:fixed))
stats2 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
  sp = RipQP.K2KrylovParams(uplo = :U, kmethod=kmethod, rhs_scale=true, #δ0 = 0.,
        form_mat = true, equilibrate = true,
          preconditioner = RipQP.LDLLowPrec(T = ftype, pos = :C, warm_start = true),
          ρ_min=ρ_min, δ_min =δ_min,
          atol_min = atol_min, rtol_min = rtol_min, atol0 = atol0, rtol0 = rtol0,
          ), 
  solve_method=RipQP.IPF(), scaling = true, history=false, presolve=false,
  ),
itol = RipQP.InputTol(max_iter=50, max_time=20.0,
ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
),
display = false,
)

header = [
  "nvar",
  "ncon",
  "nvarS",
  "pdd",
  "Pfeas",
  "Dfeas",
  "timeIPM",
  "timeF",
  "timeK",
  "it",
  "Kiter",
  "nKv",
  "nnzF",
]
nh = length(header)
pbs = ["AFIRO", "AGG", "BNL1", "CRE-A", "CYCLE", "PDS-02"]
npbs = length(pbs)
data = Matrix{Any}(undef, npbs, nh)
for i in 1:npbs
  pb = pbs[i]
  qm = QuadraticModel(readqps(string(path_pb, "\\", pb, ".SIF"), mpsformat=:fixed))
  stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                  sp = RipQP.K2KrylovParams(uplo = :U, kmethod=kmethod, rhs_scale=true, #δ0 = 0.,
                        form_mat = true, equilibrate = true,
                          preconditioner = RipQP.LDLLowPrec(T = ftype, pos = :C, warm_start = true),
                          ρ_min=ρ_min, δ_min = δ_min,
                          atol_min = atol_min, rtol_min = rtol_min, atol0 = atol0, rtol0 = rtol0,
                          ), 
                  solve_method=RipQP.IPF(), scaling = true, history=true, presolve=false,
                  ),
                itol = RipQP.InputTol(max_iter=max_iter, max_time=max_time,
                ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                ),
                display = false,
                )
  reset_timer!(RipQP.to)
  stats2 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                sp = RipQP.K2KrylovParams(uplo = :U, kmethod=kmethod, rhs_scale=true, #δ0 = 0.,
                      form_mat = true, equilibrate = true,
                        preconditioner = RipQP.LDLLowPrec(T = ftype, pos = :C, warm_start = true),
                        ρ_min=ρ_min, δ_min = δ_min,
                        atol_min = atol_min, rtol_min = rtol_min, atol0 = atol0, rtol0 = rtol0,
                        ), 
                solve_method=RipQP.IPF(), scaling = true, history=false, presolve=false,
                ),
              itol = RipQP.InputTol(max_iter=max_iter, max_time=max_time,
              ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
              ),
              display = false,
              )
  TimerOutputs.complement!(RipQP.to)
  # show(RipQP.to, sortby = :firstexec)
  dictto = TimerOutputs.todict(RipQP.to)

  data[i, :] .= Any[
    qm.meta.nvar,
    qm.meta.ncon,
    stats1.solver_specific[:nvar_slack],
    stats1.solver_specific[:pdd],
    stats1.primal_feas,
    stats1.dual_feas,
    stats2.elapsed_time,
    dictto["inner_timers"]["ripqp"]["inner_timers"]["solver IPF"]["inner_timers"]["preconditioner update"]["inner_timers"]["LDL factorize"]["time_ns"] / 1.0e9,
    dictto["inner_timers"]["ripqp"]["inner_timers"]["solver IPF"]["inner_timers"]["Krylov solve"]["time_ns"] / 1.0e9,
    stats1.solver_specific[:absolute_iter_cnt],
    sum(stats1.solver_specific[:nprodH]),
    max(20, maximum(stats1.solver_specific[:nprodH])),
    stats1.solver_specific[:nnzLDL],
    ]
    println(pb, " done")
end

pretty_table(data; 
  header = header,
  row_names= pbs,
  title = "test",
  backend = Val(:latex),
  formatters = ft_printf("%7.1e", [4, 5, 6, 7, 8, 9]))