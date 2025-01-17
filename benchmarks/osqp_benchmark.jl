include("/home/mgi.polymtl.ca/geleco/git_workspace/docGL/solvers/osqp-qm.jl")
using SolverBenchmark, JLD2

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end

function osqp2(QM)
  stats = osqp(QM, eps_abs=1.0e-3, eps_rel=1.0e-3, time_limit=1200., max_iter=80000, verbose=false)
end

path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb_lp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
# path_pb_qp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
pb2 = string(path_pb_qp, "/DUAL1.SIF")
qpdata2 = readqps(pb2);
qm2 = createQuadraticModel(qpdata2)
stats2 = osqp2(qm2)  # compile code

function optimize_osqp(path_pb)
    problems = []
    i_max = 1000
    i = 1
    for file_name in readdir(path_pb)
        if file_name[end-3:end] == ".SIF" && !(file_name in ["80BAU3B.SIF" ; "BORE3D.SIF";
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
            push!(problems, createQuadraticModel(qpdata_i, name=file_name[1:end-4]))

            if i==i_max
                break
            end
            i += 1
        end
    end
    return solve_problems(osqp2, problems)
end

save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/results"
# save_path = "/user/eleves/gleconte2017/Z/Documents/TFE/results"
# save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\results"

problems_stats_lp =  optimize_osqp(path_pb_lp)

file_lp = jldopen(string(save_path, "/osqp_lp1.jld2"), "w")
file_lp["stats"] = problems_stats_lp
close(file_lp)

problems_stats_qp =  optimize_osqp(path_pb_qp)

file_qp = jldopen(string(save_path, "/osqp_qp1.jld2"), "w")
file_qp["stats"] = problems_stats_qp
close(file_qp)
