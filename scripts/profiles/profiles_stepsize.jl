using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, BenchmarkProfiles, Plots, PrettyTables
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
# using RipQP

path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
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

function ripqpK2(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2KrylovParams(kmethod=:minres, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true#, stepsize = stepsize,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0,
                ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                ))
end

function ripqpK2_5(qm)
  return RipQP.ripqp(qm, display = false, iconf = RipQP.InputConfig(
                    sp = RipQP.K2KrylovParams(kmethod=:minres, preconditioner = :Identity, atol_min=1.0e-10, rtol_min=1.0e-10), 
                    solve_method=:IPF, history=true, stepsize = :NInf,
                    # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                    ),
                itol = RipQP.InputTol(max_iter=30, max_time=20.0,
                ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                ))
end

function optimize_ripqp!(path_pb :: String, ripqp_func :: Function, n_pb :: Int, data, row)
  problems = []
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
           push!(problems, qm)
           stats = ripqp_func(qm)
           data[row, i] = minimum(stats.solver_specific[:pddH])
           if i==i_max
               break
           end
           i += 1
       end
  end
end

n_pb = 10
n_solvers = 2
data = zeros(2, n_pb)
solvers = [:ripqpK2, :ripqpK2_5]
for solver in solvers
  optimize_ripqp!(path_pb, eval(solver), n_pb, data, 1)
  optimize_ripqp!(path_pb, eval(solver), n_pb, data, 2)
end
perf = performance_profile(PlotsBackend(), data, [string(solver) for solver in solvers], legend=:bottomright)
# plot(perf, )
title!("Performance profile (Netlib problems)")