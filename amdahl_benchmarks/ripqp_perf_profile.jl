using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, BenchmarkProfiles, Plots
using DelimitedFiles, JLD2

# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using RipQP

# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
# save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks"
save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/perf_profiles/lp2"
# save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/perf_profiles/test_qp2"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
# qm = QuadraticModel(readqps(string(path_pb, "\\irish-electricity.mps")))

function createQuadraticModel(qpdata; name="qp_pb")
  return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
          Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
          lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
          c0=qpdata.c0, name=name)
end

function zeros_logscale!(v, min_val)
  for i=1:length(v)
    if v[i] == 0
      v[i] += min_val
    end
  end
end

function K2_LDL(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
    sp = RipQP.K2LDLParams(), 
    solve_method=:IPF, #, stepsize = stepsize,
    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
    ),
  itol = RipQP.InputTol(max_iter=100, max_time=300.0,
                        ϵ_pdd = 1.0e-4, ϵ_rb = 1.0e-4, ϵ_rc = 1.0e-4,
                        ))
end

# init solvers list
fsolvers = [K2_LDL]
solvers = [:K2_LDL]

function push_solver!(
  fsolvers::Vector{Function},
  solv_str::Symbol
  formulation::Symbol,
  kmethod::Symbol,
  preconditioner::Symbol,
)
  if kmethod ∈ [:tricg, trimr, :gpmr, :lslq, :lsqr, :lsmr, :lnlq, :craig, :craigmr]
    sp = eval(formulation)(kmethod = kmethod,
                           atol_min=1.0e-10, rtol_min=1.0e-10,
                           ρ_min = 1e2 * sqrt(eps()), δ_min = 1e2 * sqrt(eps()),
                           mem = 20,
                           )
  else
    sp = eval(formulation)(kmethod = kmethod, preconditioner = preconditioner, 
                           atol_min=1.0e-10, rtol_min=1.0e-10,
                           ρ_min = 1e2 * sqrt(eps()), δ_min = 1e2 * sqrt(eps()),
                           mem = 20,
                           )
  end
  @eval function $solv_str(qm)
    return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                      sp = sp, 
                      solve_method=:IPF, #, stepsize = stepsize,
                      # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                      ),
                  itol = RipQP.InputTol(max_iter=100, max_time=300.0,
                                        ϵ_pdd = 1.0e-4, ϵ_rb = 1.0e-4, ϵ_rc = 1.0e-4,
                                        ))
  end
  push!(fsolvers, eval(solv_str))
end

# K1
for kmethod ∈ [:cg, :cg_lanczos, :cr, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K1_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K1KrylovParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K1.1 Structured
for kmethod ∈ [:lslq, :lsqr, :lsmr]
  solv_str = Symbol(:K1_1_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K1_1StructuredParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K1.2 Structured
for kmethod ∈ [:lnlq, :craig, :craigmr]
  solv_str = Symbol(:K1_2_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K1_2StructuredParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K2
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K2_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K2KrylovParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K2 Jacobi
for kmethod ∈ [:minres, :minres_qlp, :symmlq]
  solv_str = Symbol(:K2Jacobi_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K2KrylovParams, kmethod, :Jacobi)
  push!(solvers, solv_str)
end

# K2 structured
for kmethod ∈ [:tricg, :trimr, :gpmr]
  str = Symbol(:K2_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K2StructuredParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K2.5
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K2_5_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K2_5KrylovParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K2.5 Jacobi
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K2_5Jacobi_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K2_5KrylovParams, kmethod, :Jacobi)
  push!(solvers, solv_str)
end

# K2.5 Structured
for kmethod ∈ [:tricg, :trimr, :gpmr]
  solv_str = Symbol(:K2_5_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K2_5StructuredParams, kmethod, :Identity)
end

# K3
for kmethod ∈ [:bilq, :bicgstab, :usymlq, :usymqr, :qmr, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K3_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K3KrylovParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K3S
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K3S_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K3SKrylovParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K3S Structured
for kmethod ∈ [:tricg, :trimr, :gpmr]
  solv_str = Symbol(:K3S_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K3SStructuredParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K3.5
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :fom, :gmres, :dqgmres]
  solv_str = Symbol(:K3_5_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K3_5KrylovParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

# K3.5 Structured
for kmethod ∈ [:tricg, :trimr, :gpmr]
  solv_str = Symbol(:K3_5_, kmethod)
  push_solver!(fsolvers, solv_str, :RipQP.K3_5StructuredParams, kmethod, :Identity)
  push!(solvers, solv_str)
end

function get_QuadraticModels(path_pb :: String, n_pb :: Int)
  qms = []
  i_max = n_pb
  i = 1
  for file_name in readdir(path_pb)
    if file_name[end-3:end] == ".SIF"
      println(file_name)
      pb_i = string(path_pb, "/", file_name)
      if file_name in ["BLEND.SIF"; "DFL001.SIF"; "FORPLAN.SIF"; "GFRD-PNC.SIF"; "SIERRA.SIF";
                  "EXDATA.SIF"; "QFORPLAN.SIF"; "QGFRDXPN.SIF"; "VALUES.SIF"]
          qpdata_i = readqps(pb_i, mpsformat=:fixed)
      else
          qpdata_i = readqps(pb_i)
      end
      qm = createQuadraticModel(qpdata_i, name=file_name[1:end-4])
      push!(qms, qm)
      if i==i_max
        break
      end
      i += 1
    end
  end
  return qms
end

function optimize_ripqp!(qms, ripqp_func :: Function, data, solver)
  n_pb = length(qms)
  n_iter = size(data)[1]
  for i=1:n_pb
    qm = qms[i]
    stats = ripqp_func(qm)
    ldata = length(stats.solver_specific[:pddH])
    data[1:ldata, i, solver] .= stats.solver_specific[:pddH]
    if ldata < n_iter
      data[ldata+1:end, i, solver] .= stats.solver_specific[:pddH][end]
    end
  end
end

n_pb = 1000

pb_i = string(path_pb, "/", "AFIRO.SIF")
qpdata_i = readqps(pb_i)
qm = createQuadraticModel(qpdata_i, name="AFIRO")
for solver in fsolvers
  stats_compile = solver(qm)
end

n_solvers = length(solvers)

function save_problems(file_path :: String, ripqp_func :: Function, qms)
  stats = solve_problems(ripqp_func, qms)
  file = jldopen(string(file_path, "/", string(ripqp_func), ".jld2"), "w")
  file["stats"] = stats
  close(file)
  return Nothing
end

qms = get_QuadraticModels(path_pb, n_pb)
for is in 1: length(fsolvers)
  save_problems(save_path, fsolvers[is], qms)
  println(string(solvers[is]), " done")
end

solvers_list = [string(solver) for solver in solvers]

open(string(save_path, "/test2_solvs.txt"), "w") do io
  writedlm(io, solvers_list)
end;