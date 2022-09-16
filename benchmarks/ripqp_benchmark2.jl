using QuadraticModels, QPSReader
using RipQP
# using QuadraticModelsGurobi, QuadraticModelsCPLEX, QuadraticModelsXpress
# using JLD2
using CSV
using SolverBenchmark
using HSL, QDLDL
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
save_path = "/home/gelecd/code/docGL/benchmarks/ripqp_paper"
# save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
pb = string(path_pb_lp, "/AFIRO.SIF")
# pb2 = string(path_pb_qp, "/DUAL1.SIF")
qpdata = readqps(pb);
qm = createQuadraticModel(qpdata)

ripqp1(QM) = ripqp(QM, sp = K2LDLParams(),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp1(qm)
# ripqp2(QM) = ripqp(QM, sp = K2LDLParams(), kc = -1,
#                     itol = InputTol(max_iter = 800, max_time=1200.))
# stats = ripqp2(qm)
# ripqpma57(QM) = ripqp(QM,
#                     sp = K2LDLParams(fact_alg = HSLMA57Fact()),
#                     itol = InputTol(max_iter = 800, max_time=1200.))
# stats = ripqpma57(qm)
# ripqpma57_multi(QM) = ripqp(QM, mode = :multi,
#                     sp = K2LDLParams(fact_alg = HSLMA57Fact()),
#                     itol = InputTol(max_iter = 100, max_time=1200.))
# stats = ripqpma57_multi(qm)

# ripqpma57_multi2(QM) = ripqp(QM, mode = :multi,
#                     sp = K2LDLParams(fact_alg = HSLMA57Fact()),
#                     itol = InputTol(max_iter = 100, max_iter32 = 5, max_time=1200.))
# stats = ripqpma57_multi2(qm)

ripqp_ldlprecond(QM) = ripqp(QM, mode = :multi,
                    sp = K2KrylovParams(uplo = :U,
                        form_mat = true, equilibrate = false, kmethod = :gmres,
                        preconditioner = LDL(T = Float32, pos = :R, warm_start = true),
                        ρ_min=1.0e-8, δ_min = 1.0e-8,
                        mem = 10,
                        itmax = 10,
                        atol0 = 1.0e-2, rtol0 = 1.0e-2,
                        atol_min = 1.0e-8, rtol_min = 1.0e-8,
                        ),
                        sp2 = K2KrylovParams(uplo = :U,
                        form_mat = true, equilibrate = false, kmethod = :gmres,
                        preconditioner = LDL(T = Float64, pos = :R, warm_start = true),
                        ρ_min=1.0e-8, δ_min = 1.0e-8,
                        mem = 5,
                        itmax = 5,
                        atol0 = 1.0e-2, rtol0 = 1.0e-2,
                        atol_min = 1.0e-10, rtol_min = 1.0e-10,
                        ),
                    solve_method = PC(),
                    itol = InputTol(max_iter = 800, max_time=1200.,
                                    ϵ_pdd1 = 1.0e-8, ϵ_rb1 = 1.0e-6,
                                    ϵ_rc1 = 1.0e-6))

# ripqp_ldlprecond(QM) = ripqp(QM, mode = :multi, 
#                     sp = K2KrylovParams(uplo = :U,
#                         form_mat = true, equilibrate = false, kmethod = :gmres,
#                         preconditioner = LDL(T = Float32, pos = :R, warm_start = true),
#                         ρ_min=1.0e-8, δ_min = 1.0e-8,
#                         mem = 10,
#                         itmax = 10,
#                         atol0 = 1.0e-2, rtol0 = 1.0e-2,
#                         atol_min = 1.0e-8, rtol_min = 1.0e-8,
#                         ),
#                         sp2 = K2LDLParams(fact_alg = LDLFact(regul = :dynamic)),
#                     solve_method = PC(),
#                     itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp_ldlprecond(qm)
# ripqpma97(QM) = ripqp(QM,
#                     sp = K2LDLParams(fact_alg = HSLMA97Fact()),
#                     itol = InputTol(max_iter = 800, max_time=1200.))
# stats = ripqpma97(qm)
# ripqpqdldl(QM) = ripqp(QM, 
#                     sp = K2LDLParams(fact_alg = QDLDLFact()),
#                     itol = InputTol(max_iter = 800, max_time=1200.))
# stats = ripqpqdldl(qm)
# ripqpcholmod(QM) = ripqp(QM, 
#                     sp = K2LDLParams(fact_alg = CholmodFact()),
#                     itol = InputTol(max_iter = 800, max_time=1200.))
# stats = ripqpcholmod(qm)
# ripqp_nops(QM) = ripqp(QM, ps = false, itol = InputTol(max_iter = 800, max_time=1200.))
# stats = ripqp_nops(qm)
# cplex2_nops(QM) = cplex(QM, presolve=0, crossover=2, display=0, threads=1)
# cplex2(QM) = cplex(QM, crossover=2, display=0, threads=1)
# stats = cplex2_nops(qm)  # compile code
# stats = cplex2(qm)  # compile code
# gurobi2_nops(QM) = gurobi(QM, presolve=0, crossover=0, display=0, threads=1)
# gurobi2(QM) = gurobi(QM, crossover=0, display=0, threads=1)
# stats = gurobi2_nops(qm)  # compile code
# stats = gurobi2(qm)  # compile code
# xpress2_nops(QM) = xpress(QM, presolve=0, crossover=0, threads=1)
# xpress2(QM) = xpress(QM, crossover=0, threads=1)
# stats = xpress2_nops(qm)  # compile code
# stats = xpress2(qm)  # compile code
ripqp_bm_multi(QM) = ripqp(QM, mode=:multi, itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp_bm_multi(qm)
# ripqp_bm_minres(QM) = ripqp(QM, iconf = InputConfig(sp = K2_5hybridParams(preconditioner = :ActiveCHybridLDL)),
#                             itol = InputTol(max_iter=400, max_time=10.) )#,

function optimize_ripqp(path_pb :: String, ripqp_func :: Function)
    problems = []
    i_max = 1000
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

# save_problems(string(save_path, "/gurobi_nops1"), gurobi2_nops)
# save_problems(string(save_path, "/gurobi1"), gurobi2)
# save_problems(string(save_path, "/cplex_nops1"), cplex2_nops)
# save_problems(string(save_path, "/cplex1"), cplex2)
# save_problems(string(save_path, "/xpress_nops1"), xpress2_nops)
# save_problems(string(save_path, "/xpress1"), xpress2)
# save_problems(string(save_path, "/ripqp_multi1"), ripqp_bm_multi)
# save_problems(string(save_path, "/ripqp3"), ripqp1)
save_problems(string(save_path, "/ripqp_ldlprecond2"), ripqp_ldlprecond)
# save_problems(string(save_path, "/ripqp_cc1"), ripqp2)
# save_problems(string(save_path, "/ripqp_ma572"), ripqpma57)
# save_problems(string(save_path, "/ripqp_ma971"), ripqpma97)
# save_problems(string(save_path, "/ripqp_ma57_multi1"), ripqpma57_multi)
# save_problems(string(save_path, "/ripqp_ma57_multi2"), ripqpma57_multi2)
# save_problems(string(save_path, "/ripqp_ma57nosqd2"), ripqpma57_nosqd)
# save_problems(string(save_path, "/ripqp_qdldl1"), ripqpqdldl)
# save_problems(string(save_path, "/ripqp_cholmod1"), ripqpcholmod)
# save_problems(string(save_path, "/ripqp_nops1"), ripqp_nops)
# save_problems(string(save_path, "\\test"), ripqp_bm_classic)
# save_problems(string(save_path, "/ripqp_presolve_1"), ripqp_bm_presolve)
# save_problems(string(save_path, "/ripqp_mono_IPFK2_3"), ripqp_bm_classic)
# save_problems(string(save_path, "/ripqp_ccorr_1"), ripqp_bm_cc)
# save_problems(string(save_path, "/ripqp_multi_z"), ripqp_bm_multizoom)
# save_problems(string(save_path, "/ripqp_multi_K2"), ripqp_bm_multi)

# df2_lp = CSV.read(string(save_path, "\\test_lp.csv"), DataFrame)
# df2_qp = CSV.read(string(save_path, "\\test_qp.csv"), DataFrame)