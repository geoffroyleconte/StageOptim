using QuadraticModels, QPSReader
using QuadraticModelsGurobi, QuadraticModelsCPLEX, QuadraticModelsXpress
using CSV
using SolverBenchmark
using HSL, QDLDL
using RipQP

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
pb = string(path_pb_lp, "/AGG.SIF")
# pb2 = string(path_pb_qp, "/AUG2D.SIF")
qpdata = readqps(pb);
qm = createQuadraticModel(qpdata)

ripqp1(QM) = ripqp(QM, sp = K2LDLParams(),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp1(qm)
ripqp_cc(QM) = ripqp(QM, sp = K2LDLParams(), kc = -1,
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp_cc(qm)
ripqpma57(QM) = ripqp(QM,
                    sp = K2LDLParams(fact_alg = HSLMA57Fact()),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqpma57(qm)

ripqpma57_multi1(QM) = ripqp(QM, mode = :multi,
                             sp = K2LDLParams{Float32}(fact_alg = HSLMA57Fact()),
                             itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqpma57_multi1(qm)

ripqpma57_multi2(QM) = ripqp(QM, mode = :multi, early_multi_stop = false,
                    sp = K2LDLParams{Float32}(safety_dist_bnd = false,
                        fact_alg = HSLMA57Fact(), ρ_min=Float32(1.0e-7), δ_min = Float32(1.0e-7)),
                    itol = InputTol(max_iter = 800, max_time=1200.,
                                    ϵ_pdd1 = 1.0e0, ϵ_rb1 = 1.0e-2, ϵ_rc1 = 1.0e-2))
stats = ripqpma57_multi2(qm)

ripqp_ldlprecond1(QM) = ripqp(QM, mode = :multi,
                    sp = K2KrylovParams(uplo = :U,
                        form_mat = true, equilibrate = false, kmethod = :gmres,
                        preconditioner = LDL(T = Float32, warm_start = true),
                        ρ_min=1.0e-8, δ_min = 1.0e-8,
                        mem = 10, itmax = 10,
                        atol0 = 1.0e-2, rtol0 = 1.0e-2,
                        atol_min = 1.0e-8, rtol_min = 1.0e-8,
                        ),
                        sp2 = K2KrylovParams(uplo = :U,
                        form_mat = true, equilibrate = false, kmethod = :gmres,
                        preconditioner = LDL(T = Float64, warm_start = true),
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
stats = ripqp_ldlprecond1(qm)
ripqp_ldlprecond2(QM) = ripqp(QM, mode = :multi,
                    sp = K2KrylovParams(uplo = :U,
                        form_mat = true, equilibrate = true, kmethod = :gmres,
                        preconditioner = LDL(T = Float32, warm_start = true),
                        ρ_min=1.0e-8, δ_min = 1.0e-8,
                        mem = 10, itmax = 10,
                        atol0 = 1.0e-2, rtol0 = 1.0e-2,
                        atol_min = 1.0e-8, rtol_min = 1.0e-8,
                        ),
                        sp2 = K2KrylovParams(uplo = :U,
                        form_mat = true, equilibrate = true, kmethod = :gmres,
                        preconditioner = LDL(T = Float64, warm_start = true),
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
stats = ripqp_ldlprecond2(qm)

ripqp_lldlprecond(QM) = ripqp(QM, mode = :multi, 
                    sp = K2KrylovParams(uplo = :L,
                        form_mat = true, equilibrate = false, kmethod = :gmres,
                        preconditioner = LDL(fact_alg = LLDLFact(mem=20), pos = :R,
                                            T = Float32, warm_start = true),
                        ρ_min=1.0e-8, δ_min = 1.0e-8,
                        mem = 10, itmax = 10,
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
                    solve_method = IPF(), solve_method2 = PC(),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp_lldlprecond(qm)


ripqpma97(QM) = ripqp(QM,
                    sp = K2LDLParams(fact_alg = HSLMA97Fact()),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqpma97(qm)
ripqpqdldl(QM) = ripqp(QM, 
                    sp = K2LDLParams(fact_alg = QDLDLFact()),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqpqdldl(qm)
ripqpcholmod(QM) = ripqp(QM, 
                    sp = K2LDLParams(fact_alg = CholmodFact()),
                    itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqpcholmod(qm)
cplex2(QM) = cplex(QM, crossover=2, display=0, threads=1)
stats = cplex2(qm)  # compile code
gurobi2(QM) = gurobi(QM, crossover=0, display=0, threads=1)
stats = gurobi2(qm)  # compile code
xpress2(QM) = xpress(QM, crossover=0, threads=1)
stats = xpress2(qm)  # compile code
ripqp_bm_multi(QM) = ripqp(QM, mode=:multi, itol = InputTol(max_iter = 800, max_time=1200.))
stats = ripqp_bm_multi(qm)

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

save_problems(string(save_path, "/ripqp1"), ripqp1)
save_problems(string(save_path, "/gurobi1"), gurobi2)
save_problems(string(save_path, "/cplex1"), cplex2)
save_problems(string(save_path, "/xpress1"), xpress2)

save_problems(string(save_path, "/ripqp_cc1"), ripqp_cc)

save_problems(string(save_path, "/ripqp_ma57"), ripqpma57)
save_problems(string(save_path, "/ripqp_ma971"), ripqpma97)
save_problems(string(save_path, "/ripqp_qdldl1"), ripqpqdldl)
save_problems(string(save_path, "/ripqp_cholmod1"), ripqpcholmod)

save_problems(string(save_path, "/ripqp_multi1"), ripqp_bm_multi)
save_problems(string(save_path, "/ripqp_ma57_multi1"), ripqpma57_multi1) # check regu
save_problems(string(save_path, "/ripqp_ma57_multi2"), ripqpma57_multi2) # check regu

save_problems(string(save_path, "/ripqp_ldlprecond1"), ripqp_ldlprecond1)
save_problems(string(save_path, "/ripqp_ldlprecond2"), ripqp_ldlprecond2)
save_problems(string(save_path, "/ripqp_lldlprecond"), ripqp_lldlprecond)

# save_problems(string(save_path, "/ripqp_ldlprecondma57"), ripqp_ldlprecond)

using Quadmath, DoubleFloats
# path_pb = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\quadLP\data\MPS"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"
path_pb = "/home/gelecd/datasets/quad_problems"
save_path = "/home/gelecd/code/docGL/benchmarks/ripqp_paper"

function createQuadraticModel_T(qpdata; T = Float128, name="qp_pb")
    return QuadraticModel(convert(Array{T}, qpdata.c), qpdata.qrows, qpdata.qcols,
            convert(Array{T}, qpdata.qvals),
            Arows=qpdata.arows, Acols=qpdata.acols,
            Avals=convert(Array{T}, qpdata.avals),
            lcon=convert(Array{T}, qpdata.lcon),
            ucon=convert(Array{T}, qpdata.ucon),
            lvar=convert(Array{T}, qpdata.lvar),
            uvar=convert(Array{T}, qpdata.uvar),
            c0=T(qpdata.c0), x0 = zeros(T, length(qpdata.c)), name=name)
end

T = Float128
Tlow = Float64

# compile
path_pb_lp = "/home/gelecd/.julia/artifacts/545f8c5577a056981a21caf3f53bd7b59cf67410/optrove-netlib-lp-f83996fca937"
qm1 = createQuadraticModel_T(readqps(string(path_pb_lp, "/AFIRO.SIF")), T=T)

ripqp_multi_quad(qm; T = T, Tlow = Tlow) = ripqp(qm, 
  mode = :multi,
  early_multi_stop = false,
  sp = K2KrylovParams{Tlow}( # solve in Float64
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 50,
    mem = 50,
    preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
    ρ_min=1.0e-15,
    δ_min = 1.0e-15,
    atol_min = 1.0e-16,
    rtol_min = 1.0e-16,
  ),
    sp2 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 5,
    mem = 5,
    preconditioner = LDL(T = T, pos = :R, warm_start = true),
    ρ_min=T(1.0e-15),
    δ_min = T(1.0e-15),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(),
  solve_method2=PC(),
  itol = InputTol(T, max_iter = 7000, max_time = 20000.0, max_iter1 = 100, ϵ_pdd1 = T(1.0e1),
    ϵ_rc1 = T(1.0e-6), ϵ_rb1 = T(1.0e-6)),
  display = true,
)
stats = ripqp_multi_quad(qm1)


function optimize_ripqp(path_pb :: String, ripqp_func :: Function, T::DataType)
  problems = [
    createQuadraticModel_T(readqps(string(path_pb, "/TMA_ME_presolved.mps")), T = T, name = "TMA_ME"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAlift_presolved.mps")), T = T, name = "GlcAlift"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAerWT_presolved.mps")), T = T, name = "GlcAerWT"),
  ]

  return solve_problems(ripqp_func, problems)
end

function save_quad_problems(file_path :: String, ripqp_func :: Function; path_pb :: String = path_pb, T = T)
  lp_stats = optimize_ripqp(path_pb, ripqp_func, T)
  CSV.write(string(file_path, "_quad.csv"), lp_stats)
  return Nothing
end

function optimize_ripqp_nops(path_pb :: String, ripqp_func :: Function, T::DataType)
  problems = [
    createQuadraticModel_T(readqps(string(path_pb, "/TMA_ME.mps")), T = T, name = "TMA_ME"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAlift.mps")), T = T, name = "GlcAlift"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAerWT.mps")), T = T, name = "GlcAerWT"),
  ]

  return solve_problems(ripqp_func, problems)
end

function save_quad_problems_nops(file_path :: String, ripqp_func :: Function; path_pb :: String = path_pb, T = T)
  lp_stats = optimize_ripqp_nops(path_pb, ripqp_func, T)
  CSV.write(string(file_path, "_nops_quad.csv"), lp_stats)
  return Nothing
end

save_quad_problems(string(save_path, "/ripqp_multi1"), ripqp_multi, T = T)
save_quad_problems(string(save_path, "/ripqp_mono1"), ripqp_mono, T = T)
save_quad_problems(string(save_path, "/ripqp_multi"), ripqp_multi_quad, T = T)
save_quad_problems_nops(string(save_path, "/ripqp_multi"), ripqp_multi_quad, T = T)