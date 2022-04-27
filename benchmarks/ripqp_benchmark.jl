using QuadraticModels, QPSReader
using RipQP
# using JLD2
using CSV
using SolverBenchmark
# include("/home/mgi.polymtl.ca/geleco/git_workspace/docGL/utils/K1QR.jl")

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end

# path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
path_pb_lp = "/home/gelecd/.julia/artifacts/545f8c5577a056981a21caf3f53bd7b59cf67410/optrove-netlib-lp-f83996fca937"
# path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
path_pb_qp = "/home/gelecd/.julia/artifacts/0eff5ae5b345db85386f55f672a19c90f23257b2/optrove-maros-meszaros-9adfb5707b1e"
# path_pb_lp = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb_qp = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
# save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/results"
save_path = "/home/gelecd/code/docGL/benchmarks/frontal22_results/bm1"
# save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
pb = string(path_pb_lp, "/AFIRO.SIF")
# pb2 = string(path_pb_qp, "/DUAL1.SIF")
qpdata = readqps(pb);
qm = createQuadraticModel(qpdata)
stats =  ripqp(qm)  # compile code

ripqp_bm_classic(QM) = ripqp(QM, itol = InputTol(max_time=200.))
# ripqp_bm_cc(QM) = ripqp(QM, iconf = InputConfig(kc=-1), itol = InputTol(max_time=1200.))
# ripqp_bm_presolve(QM) =  ripqp(QM, itol = InputTol(max_time=1200.), iconf = InputConfig(presolve=true, scaling=true))
# ripqp_bm_multiref(QM) = ripqp(QM, iconf = InputConfig(mode=:multi, refinement=:multiref), itol = InputTol(max_time=1200.))
# ripqp_bm_multi(QM) = ripqp(QM, iconf = InputConfig(mode=:multi), itol = InputTol(max_time=1200.))
# ripqp_bm_minres(QM) = ripqp(QM, iconf = InputConfig(sp = K2_5hybridParams(preconditioner = :ActiveCHybridLDL)),
#                             itol = InputTol(max_iter=400, max_time=10.) )#,

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

function save_problems(file_path :: String, ripqp_func :: Function, 
                       path_pb_lp :: String = path_pb_lp, path_pb_qp :: String = path_pb_qp)

    lp_classic =  optimize_ripqp(path_pb_lp, ripqp_func)
    CSV.write(string(file_path, "_lp.csv"), lp_classic)
    qp_classic =  optimize_ripqp(path_pb_qp, ripqp_func)
    CSV.write(string(file_path, "_qp.csv"), qp_classic)
    
    return Nothing
end

save_problems(string(save_path, "/ripqp_mono_wtimerout"), ripqp_bm_classic)
# save_problems(string(save_path, "\\test"), ripqp_bm_classic)
# save_problems(string(save_path, "/ripqp_presolve_1"), ripqp_bm_presolve)
# save_problems(string(save_path, "/ripqp_mono_IPFK2_3"), ripqp_bm_classic)
# save_problems(string(save_path, "/ripqp_ccorr_1"), ripqp_bm_cc)
# save_problems(string(save_path, "/ripqp_multi"), ripqp_bm_multi)
# save_problems(string(save_path, "/ripqp_multi_r"), ripqp_bm_multiref)
# save_problems(string(save_path, "/ripqp_mono_z"), ripqp_bm_zoom)
# save_problems(string(save_path, "/ripqp_multi_z"), ripqp_bm_multizoom)
# save_problems(string(save_path, "/ripqp_multi_K2"), ripqp_bm_multi)

# df2_lp = CSV.read(string(save_path, "\\test_lp.csv"), DataFrame)
# df2_qp = CSV.read(string(save_path, "\\test_qp.csv"), DataFrame)