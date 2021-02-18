using QuadraticModels, QPSReader
using RipQP
using JLD2
using SolverBenchmark

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end

path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb_lp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
# path_pb_qp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
pb = string(path_pb_lp, "/AFIRO.SIF")
# pb2 = string(path_pb_qp, "/DUAL1.SIF")
qpdata = readqps(pb);
qm = createQuadraticModel(qpdata)
stats =  ripqp(qm, mode=:mono, regul=:classic, K=0)  # compile code

ripqp_bm_classic(QM) = ripqp(QM, mode=:mono, regul=:classic, K=0, max_time=1.)
ripqp_bm_classic_nrtol(QM) = ripqp(QM, mode=:mono, regul=:classic, K=0, max_time=1. , normalize_rtol=false)
ripqp_bm_dynamic(QM) = ripqp(QM, mode=:mono, regul=:dynamic, K=0, max_time=1.)
ripqp_bm_ccorr(QM) = ripqp(QM, mode=:mono, regul=:classic, K=-1, max_time=1.)
ripqp_bm_multi(QM) = ripqp(QM, mode=:multi, regul=:classic, K=0, max_time=1.)
ripqp_bm_muti_dynamic(QM) = ripqp(QM, mode=:multi, regul=:classic, K=0, max_time=1.)

function optimize_ripqp(path_pb :: String, ripqp_func :: Function)
    problems = []
    i_max = 1000
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
             push!(problems, createQuadraticModel(qpdata_i, name=file_name[1:end-4]))

             if i==i_max
                 break
             end
             i += 1
         end
    end

    return solve_problems(ripqp_func, problems)
end


save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/StageOptim/amdahl_benchmarks/results"
# save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\StageOptim\\amdahl_benchmarks\\results"

function save_problems(file_path :: String, ripqp_func :: Function, 
                       path_pb_lp :: String = path_pb_lp, path_pb_qp :: String = path_pb_qp)

    lp_classic =  optimize_ripqp(path_pb_lp, ripqp_func)
    file_lp = jldopen(string(file_path, "_lp.jld2"), "w")
    file_lp["stats"] = lp_classic
    close(file_lp)
    qp_classic = optimize_ripqp(path_pb_qp, ripqp_func)
    file_qp = jldopen(string(file_path, "_qp.jld2"), "w")
    file_qp["stats"] = qp_classic
    close(file_qp)
    
    return Nothing
end

save_problems(string(save_path, "/ripqp_mono_1"), ripqp_bm_classic)
save_problems(string(save_path, "/ripqp_mono_nrtol_1"), ripqp_bm_classic_nrtol)
save_problems(string(save_path, "/ripqp_dynamic_1"), ripqp_bm_dynamic)
save_problems(string(save_path, "/ripqp_ccorr_1"), ripqp_bm_ccorr)
save_problems(string(save_path, "/ripqp_multi_1"), ripqp_bm_multi)
save_problems(string(save_path, "/ripqp_multi_dynamic_1"), ripqp_bm_muti_dynamic)
