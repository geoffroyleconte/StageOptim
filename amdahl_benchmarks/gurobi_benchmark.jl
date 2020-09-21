using QuadraticModelsGurobi
using QPSReader
using QuadraticModels
using NLPModels
using SolverTools
using SolverBenchmark
using DataFrames
using LinearAlgebra
using JLD2
using SparseArrays

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end

function gurobi2(QM)
    return gurobi(QM, presolve=0, crossover=0, display=0)
end

path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb = "/user/eleves/gleconte2017/Z/Documents/TFE/marosmeszaros"
#path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
pb2 = string(path_pb_qp, "/DUAL1.SIF")
qpdata2 = readqps(pb2);
qm2 = createQuadraticModel(qpdata2)
stats2 = gurobi(qm2)  # compile code

function optimize_gurobi(path_pb)
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
    return solve_problems(gurobi2, problems)
end

save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/StageOptim/amdahl_benchmarks/results"
# save_path = "/user/eleves/gleconte2017/Z/Documents/TFE/results"
#save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\results"

problems_stats_lp =  optimize_gurobi(path_pb_lp)

file_lp = jldopen(string(save_path, "/gurobi_scaling_lp3.jld2"), "w")
file_lp["stats"] = problems_stats_lp
close(file_lp)

problems_stats_qp =  optimize_gurobi(path_pb_qp)

file_qp = jldopen(string(save_path, "/gurobi_scaling_qp3.jld2"), "w")
file_qp["stats"] = problems_stats_qp
close(file_qp)
