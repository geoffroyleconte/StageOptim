# res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\ripqp_paper"
using Plots
using DataFrames, SolverBenchmark, SolverTools
using JLD2
using CSV
# using FileIO
open_file(fname; res_path = res_path) = CSV.read(string(res_path, "\\", fname, ".csv"), DataFrame)

gurobi1_lp = open_file("gurobi1_lp")
cplex1_lp = open_file("cplex1_lp")
xpress1_lp = open_file("xpress1_lp")
gurobi_nops1_lp = open_file("gurobi_nops1_lp")
cplex_nops1_lp = open_file("cplex_nops1_lp")
xpress_nops1_lp = open_file("xpress_nops1_lp")
ripqp1_lp = open_file("ripqp1_lp") # compare commercial solvs + cc1
ripqp2_lp = open_file("ripqp2_lp") # compare factos
ripqp3_lp = open_file("ripqp3_lp")  # compare multi
ripqp_ma57_lp = open_file("ripqp_ma571_lp")
ripqp_ma572_lp = open_file("ripqp_ma571_lp") # compare multi
ripqp_ma97_lp = open_file("ripqp_ma971_lp")
ripqp_ma57_multi_lp = open_file("ripqp_ma57_multi1_lp")
ripqp_ma57_multi2_lp = open_file("ripqp_ma57_multi2_lp") # 5 itmax 32
ripqp_ma57nosqd_lp = open_file("ripqp_ma57nosqd1_lp")
ripqp_ma57_lp2 = open_file("ripqp_ma572_lp") # regu comm dhab
ripqp_ma57nosqd_lp2 = open_file("ripqp_ma57nosqd2_lp")
ripqp_qdldl_lp = open_file("ripqp_qdldl1_lp")
ripqp_cholmod_lp = open_file("ripqp_cholmod1_lp")
ripqp_multi1_lp = open_file("ripqp_multi1_lp")
ripqp_multi2_lp = open_file("ripqp_multi2_lp")
ripqp_nops1_lp = open_file("ripqp_nops1_lp")
ripqp_cc1_lp = open_file("ripqp_cc1_lp")
ripqp_ldlprecond1_lp = open_file("ripqp_ldlprecond1_lp") # regu1 1.0e-8, stop crit 64, no equi
ripqp_ldlprecond2_lp = open_file("ripqp_ldlprecond2_lp") # regu1 1.0e-8 equi
ripqp_ldlprecond3_lp = open_file("ripqp_ldlprecond3_lp") # regu1 1.0e-8 stop crit 64, equi
ripqp_ldlprecond4_lp = open_file("ripqp_ldlprecond4_lp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch
ripqp_ldlprecond2_5_lp = open_file("ripqp_ldlprecond2_5_lp") 
ripqp_lldlprecond1_lp = open_file("ripqp_ldlprecond1_lp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch
ripqp_lldlprecond2_lp = open_file("ripqp_ldlprecond2_lp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch

easy_pbs_lp = findall(ripqp1_lp.elapsed_time .≤ 200.0)
stats_lp = Dict(
                # :gurobi => gurobi1_lp,
                # :cplex => cplex1_lp,
                # :xpress => xpress1_lp,
                # :gurobi_nops1 => gurobi_nops1_lp,
                # :cplex_nops1 => cplex_nops1_lp,
                # :xpress_nops1 => xpress_nops1_lp,
                # :ripqp => ripqp1_lp,
                # :ripqp_ldl => ripqp2_lp,
                # :ripqp_qdldl => ripqp_qdldl_lp,
                # :ripqp => ripqp3_lp,
                # :ripqp_ldlfact => ripqp1_lp,
                # :ripqp_multi => ripqp_multi1_lp,
                # :ripqp_multi2 => ripqp_multi2_lp,
                # :ripqp_nops1 => ripqp_nops1_lp,
                #  :ripqp_ma57 => ripqp_ma57_lp,
                #  :ripqp_ma97 => ripqp_ma97_lp,
                # :ripqp_ma57 => ripqp_ma57_lp2,
                # :ripqp_ma57_multi => ripqp_ma57_multi_lp,
                # :ripqp_ma57_multi2 => ripqp_ma57_multi2_lp,
                # :ripqp => filter(x -> x.id ∉ easy_pbs_lp, ripqp1_lp),
                # :ripqp_ma97 => filter(x -> x.id ∉ easy_pbs_lp, ripqp_ma97_lp),
                # :ripqp_ma57 => filter(x -> x.id ∉ easy_pbs_lp, ripqp_ma572_lp),
                # :ripqp_ma57_multi => filter(x -> x.id ∉ easy_pbs_lp, ripqp_ma57_multi_lp),
                # :ripqp_ma57nosqd => ripqp_ma57nosqd_lp2,
                # :ripqp_cholmod => ripqp_cholmod_lp,
                # :ripqp_cc => ripqp_cc1_lp,
                # :ripqp_multifact1 => ripqp_ldlprecond1_lp,
                # :ripqp_ldlprecond1 => ripqp_ldlprecond3_lp,
                :ripqp_ldlprecond4 => ripqp_ldlprecond4_lp,
                # :ripqp_multifact2 => ripqp_ldlprecond3_lp,
                # :ripqp_ldlprecond2_5 => ripqp_ldlprecond2_5_lp,
                )

gurobi1_qp = open_file("gurobi1_qp")
cplex1_qp = open_file("cplex1_qp")
xpress1_qp = open_file("xpress1_qp") 
gurobi_nops1_qp = open_file("gurobi_nops1_qp") 
cplex_nops1_qp = open_file("cplex_nops1_qp")
xpress_nops1_qp = open_file("xpress_nops1_qp") 
ripqp1_qp = open_file("ripqp1_qp")
ripqp2_qp = open_file("ripqp2_qp")
ripqp3_qp = open_file("ripqp3_qp")
ripqp_ma57_qp = open_file("ripqp_ma571_qp")
ripqp_ma97_qp = open_file("ripqp_ma971_qp")
ripqp_ma57_multi_qp = open_file("ripqp_ma57_multi1_qp")
ripqp_ma57_multi2_qp = open_file("ripqp_ma57_multi2_qp") # 5 itmax 32
ripqp_ma57nosqd_qp = open_file("ripqp_ma57nosqd1_qp")
ripqp_ma572_qp = open_file("ripqp_ma572_qp")
ripqp_ma57nosqd_qp2 = open_file("ripqp_ma57nosqd2_qp") 
ripqp_qdldl_qp = open_file("ripqp_qdldl1_qp") 
ripqp_cholmod_qp = open_file("ripqp_cholmod1_qp") 
ripqp_multi1_qp = open_file("ripqp_multi1_qp") 
ripqp_multi2_qp = open_file("ripqp_multi2_qp")
ripqp_nops1_qp = open_file("ripqp_nops1_qp") 
ripqp_cc1_qp = open_file("ripqp_cc1_qp")
ripqp_ldlprecond1_qp = open_file("ripqp_ldlprecond1_qp") # regu1 1.0e-8, stop crit 64, no equi
ripqp_ldlprecond2_qp = open_file("ripqp_ldlprecond2_qp") # regu1 1.0e-8 equi
ripqp_ldlprecond3_qp = open_file("ripqp_ldlprecond3_qp") # regu1 1.0e-8 stop crit 64, equi
ripqp_ldlprecond4_qp = open_file("ripqp_ldlprecond4_qp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch
ripqp_ldlprecond2_5_qp = open_file("ripqp_ldlprecond2_5_qp") 
ripqp_lldlprecond1_qp = open_file("ripqp_ldlprecond1_qp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch
ripqp_lldlprecond2_qp = open_file("ripqp_ldlprecond2_qp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch

easy_pbs_qp = findall(ripqp_ma57_qp.elapsed_time .≤ 5.0)
stats_qp = Dict(
                # :gurobi => gurobi1_qp,
                # :cplex => cplex1_qp,
                # :xpress => xpress1_qp,
                # :gurobi_nops1 => gurobi_nops1_qp,
                # :cplex_nops1 => cplex_nops1_qp,
                # :xpress_nops1 => xpress_nops1_qp,
                # :ripqp => ripqp1_qp,
                # :ripqp_ldl => ripqp2_qp,
                # :ripqp => ripqp3_qp,
                # :ripqp_ma57 => ripqp_ma572_qp,
                # :ripqp_qdldl => ripqp_qdldl_qp,
                # :ripqp_ma57 => ripqp_ma57_qp,
                # :ripqp_ma97 => ripqp_ma97_qp,
                # :ripqp_ma57 => ripqp_ma572_qp,
                # :ripqp_ma57_multi => ripqp_ma57_multi_qp,
                # :ripqp_multi2 => ripqp_multi2_qp,
                # :ripqp_ma57_multi2 => ripqp_ma57_multi2_qp,
                # :ripqp_ma57nosqd => ripqp_ma57nosqd_qp,
                # :ripqp => filter(x -> x.id ∉ easy_pbs_qp, ripqp_ma572_qp),
                # :ripqp_qdldl => filter(x -> x.id ∉ easy_pbs_qp, ripqp_qdldl_qp),
                # :ripqp_ma57_multi => filter(x -> x.id ∉ easy_pbs_qp, ripqp_ma57_multi_qp),
                # :ripqp_cholmod => ripqp_cholmod_qp,
                # :ripqp_multi => ripqp_multi1_qp,
                # :ripqp_nops1 => ripqp_nops1_qp,
                # :ripqp_cc => ripqp_cc1_qp,
                :ripqp_multifact1 => ripqp_ldlprecond1_qp,
                # :ripqp_ldlprecond2 => ripqp_ldlprecond2_qp,
                # :ripqp_multifact2 => ripqp_ldlprecond3_qp,
                :ripqp_multifact4 => ripqp_ldlprecond4_qp,
                )

function dfstat(df)
  output = zeros(length(df.status))
  for i=1:length(df.status)
    if df.primal_feas[i] === missing || df.objective[i] == Inf
      output[i] = Inf
    else 
      # output[i] = df.iter[i]
      output[i] = df.relative_iter_cnt[i]
      # output[i] = df.iters_sp[i]
      # output[i] = df.elapsed_time[i]
      # output[i] = 4 * df.iters_sp2[i] + df.iters_sp[i]
    end
    if df.status[i] ∉ ["first_order", "acceptable"]
      output[i] = Inf
    end
  end
  return output
end

cost = df -> df.elapsed_time + (df.status .!= :first_order) * Inf # + (df.elapsed_time .>= 10.) * Inf
pgfplotsx()
perf = performance_profile(stats_lp, dfstat,legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# title!("Performance profile (Netlib problems)")
# perf = performance_profile(stats_qp, dfstat,legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# title!("Performance profile (Maros and Meszaros problems)")
# display("image/svg+xml", perf)
# savefig(perf, raw"C:\Users\Geoffroy Leconte\Documents\doctorat\biblio\papiers\ripqp\paper\profiles\multifact_mm_riter.tikz")