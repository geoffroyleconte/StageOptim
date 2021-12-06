using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, BenchmarkProfiles, Plots
using DelimitedFiles, JLD2

# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using RipQP

# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
# save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks"
save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/perf_profiles/test1"
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

function ripqpLDL(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
    sp = RipQP.K2LDLParams(), 
    solve_method=:IPF, history=true#, stepsize = stepsize,
    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
    ),
  itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK1(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K1KrylovParams(kmethod=:cg, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true#, stepsize = stepsize,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK2(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2KrylovParams(kmethod=:minres, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true#, stepsize = stepsize,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK2_5(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2_5KrylovParams(kmethod=:minres, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK3(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K3KrylovParams(kmethod=:bicgstab, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK3_5(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K3_5KrylovParams(kmethod=:minres, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpTricg(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2StructuredParams(kmethod=:tricg, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpTricgK2_5(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2_5StructuredParams(kmethod=:tricg, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpTrimr(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2StructuredParams(kmethod=:trimr, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpTrimrK2_5(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2_5StructuredParams(kmethod=:trimr, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK2Jacobi(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2KrylovParams(kmethod=:minres, preconditioner = :Jacobi, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true#, stepsize = stepsize,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function ripqpK2_5Jacobi(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2_5KrylovParams(kmethod=:minres, preconditioner = :Jacobi, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0))
end

function get_QuadraticModels(path_pb :: String, n_pb :: Int)
  qms = []
  i_max = n_pb
  i = 1
  for file_name in readdir(path_pb)
    if file_name[end-3:end] == ".SIF" && !(file_name in["80BAU3B.SIF" ; "BORE3D.SIF";
                                                    "CAPRI.SIF"; "CZPROB.SIF";
                                                    "ETAMACRO.SIF"; "FINNIS.SIF";
                                                    "FORPLAN.SIF"; "GREENBEA.SIF";
                                                    "GREENBEB.SIF"; "MAROS.SIF";
                                                    "NESM.SIF"; "PEROLD.SIF";
                                                    "PILOT-JA.SIF"; "PILOT-WE.SIF";
                                                    "PILOT.SIF"; "PILOT4.SIF";
                                                    "PILOT87.SIF"; "PILOTNOV.SIF";
                                                      "RECIPELP.SIF"; "SHELL.SIF";
                                                    "SIERRA.SIF"; "STAIR.SIF";
                                                    "STANDATA.SIF"; "STANDGUB.SIF";
                                                    "STANDMPS.SIF"; "TUFF.SIF";
                                                    "VTP-BASE.SIF"; "DTOC3.SIF";
                                                    "HS35MOD.SIF";"QBORE3D.SIF";
                                                    "QCAPRI.SIF"; "QETAMACR.SIF";
                                                      "QFORPLAN.SIF"; "QPCSTAIR.SIF";
                                                    "QPCSTAIR.SIF"; "QPILOTNO.SIF";
                                                    "QRECIPE.SIF"; "QSHELL.SIF";
                                                    "QSIERRA.SIF"; "QSTAIR.SIF";
                                                    "QSTANDAT.SIF"; "UBH1.SIF";
                                                    "YAO.SIF"]) # problems with fixed variables


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

n_k = 31
n_pb = 1
solvers = [
  :ripqpLDL,
  :ripqpK1,
  :ripqpK2,
  # :ripqpK2Jacobi,
  # :ripqpK2_5,
  # :ripqpK2_5Jacobi,
  # :ripqpK3,
  # :ripqpK3_5,
  # :ripqpTricg,
  # :ripqpTricgK2_5,
  # :ripqpTrimr,
  # :ripqpTrimrK2_5,
  ]

pb_i = string(path_pb, "/", "AFIRO.SIF")
qpdata_i = readqps(pb_i)
qm = createQuadraticModel(qpdata_i, name="AFIRO")
for solver in solvers
  stats_compile = eval(solver)(qm)
end

n_solvers = length(solvers)

function save_problems(file_path :: String, ripqp_func :: Function, qms, path_pb :: String)
  stats = solve_problems(ripqp_func, qms)
  file = jldopen(string(file_path, "/", string(ripqp_func), ".jld2"), "w")
  file["stats"] = stats
  close(file)
  return Nothing
end

qms = get_QuadraticModels(path_pb, n_pb)
for is in 1: length(solvers)
  save_problems(save_path, eval(solvers[is]), qms, path_pb)
  println(string(solvers[is]), " done")
end

solvers_list = [string(solver) for solver in solvers]

open(string(save_path, "/test2_solvs.txt"), "w") do io
  writedlm(io, solvers_list)
end;