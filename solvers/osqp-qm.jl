using OSQP, QuadraticModels, SolverCore, SparseArrays, LinearAlgebra, QPSReader

const status_map_qm = Dict{Symbol,Symbol}(
  :Dual_infeasible_inaccurate => :infeasible,
  :Primal_infeasible_inaccurate => :infeasible,
  :Solved_inaccurate => :unknown,
  :Solved => :acceptable,
  :Max_iter_reached => :max_iter,
  :Primal_infeasible => :infeasible,
  :Dual_infeasible => :infeasible,
  :Interrupted => :user,
  :Time_limit_reached => :max_time,
  :Non_convex => :exception,
  :Unsolved => :unknown)

function osqp(QM::QuadraticModel; kwargs...)

  m = OSQP.Model()
  
  H = sparse(QM.data.Hcols, QM.data.Hrows, QM.data.Hvals, QM.meta.nvar, QM.meta.nvar)
  A = [sparse(QM.data.Arows, QM.data.Acols, QM.data.Avals, QM.meta.ncon, QM.meta.nvar); I]
  l = [QM.meta.lcon; QM.meta.lvar]
  u = [QM.meta.ucon; QM.meta.uvar]

  OSQP.setup!(m; P=H, q=QM.data.c, A=A, l=l, u=u, kwargs...)
  results = OSQP.solve!(m)

  stats = GenericExecutionStats(
    status_map_qm[results.info.status],
    QM,
    solution = results.x,
    objective = results.info.obj_val + QM.data.c0,
    dual_feas = results.info.dua_res,
    primal_feas = results.info.pri_res,
    multipliers = results.y,
    iter = results.info.iter,
    elapsed_time = results.info.run_time,
    solver_specific = Dict(:setup_time => results.info.setup_time,
                           :solve_time => results.info.solve_time,
                           :supdate_time => results.info.update_time,
                           :polish_time => results.info.polish_time,
                           ),
  )
  return stats
end

# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# qm = QuadraticModel(readqps(string(path_pb, "\\LISWET1.SIF"), mpsformat=:fixed))
# stats = osqp(qm, eps_abs=0.0, eps_rel=1.0e-7, time_limit=10., max_iter=8000)
# println(stats)