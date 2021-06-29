using QPALM, QuadraticModels, SolverCore, SparseArrays, LinearAlgebra, QPSReader

const status_map_qm = Dict{Symbol, Symbol}(
     :Error => :exception,
     :Solved => :acceptable,
     :Dual_terminated => :unknown,
     :Max_iter_reached => :max_iter,
     :Primal_infeasible => :infeasible,
     :Dual_infeasible => :infeasible,
     :Time_limit_reached => :max_time,
     :Unsolved => :unknown,
)

function createQuadraticModel(qpdata; name="qp_pb")
  return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
          Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
          lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
          c0=qpdata.c0, name=name)
end

function qpalm(QM::QuadraticModel; settings...)
  H = sparse(QM.data.Hcols, QM.data.Hrows, QM.data.Hvals, QM.meta.nvar, QM.meta.nvar)
  A = [sparse(QM.data.Arows, QM.data.Acols, QM.data.Avals, QM.meta.ncon, QM.meta.nvar); I]
  l = [QM.meta.lcon; QM.meta.lvar]
  u = [QM.meta.ucon; QM.meta.uvar]
  model = QPALM.Model()
  QPALM.setup!(model, Q=H, q=QM.data.c, A=A, bmin=l, bmax=u; settings...)
  results = QPALM.solve!(model)
  stats = GenericExecutionStats(
    status_map_qm[results.info.status],
    QM,
    solution = results.x,
    objective = results.info.objective + QM.data.c0,
    dual_feas = results.info.dua_res_norm,
    primal_feas = results.info.pri_res_norm,
    multipliers = results.y,
    iter = results.info.iter,
    elapsed_time = results.info.run_time,
    solver_specific = Dict(:setup_time => results.info.setup_time,
                        :solve_time => results.info.solve_time,
                        ),
  )
  return stats
end

path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb_lp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
# path_pb_qp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
pb2 = string(path_pb_qp, "/QSEBA.SIF")
qpdata2 = readqps(pb2);
qm2 = createQuadraticModel(qpdata2)
stats2 = qpalm(qm2, eps_rel=1.0e-6,eps_abs=0.0, verbose=0, time_limit=1200.) 
println(stats2)