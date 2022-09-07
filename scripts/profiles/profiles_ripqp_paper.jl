# res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\ripqp_paper"
using Plots
using DataFrames, SolverBenchmark, SolverTools
using JLD2
using CSV
# using FileIO
function open_file(path)
    file = jldopen(path, "r")
    df_out = file["stats"]
    close(file)
    return df_out
end

gurobi1_lp = CSV.read(string(res_path, "\\gurobi1_lp.csv"), DataFrame)
cplex1_lp = CSV.read(string(res_path, "\\cplex1_lp.csv"), DataFrame)
xpress1_lp = CSV.read(string(res_path, "\\xpress1_lp.csv"), DataFrame)
gurobi_nops1_lp = CSV.read(string(res_path, "\\gurobi_nops1_lp.csv"), DataFrame)
cplex_nops1_lp = CSV.read(string(res_path, "\\cplex_nops1_lp.csv"), DataFrame)
xpress_nops1_lp = CSV.read(string(res_path, "\\xpress_nops1_lp.csv"), DataFrame)
ripqp1_lp = CSV.read(string(res_path, "\\ripqp1_lp.csv"), DataFrame) # compare commercial solvs + cc1
ripqp2_lp = CSV.read(string(res_path, "\\ripqp2_lp.csv"), DataFrame) # compare factos
ripqp3_lp = CSV.read(string(res_path, "\\ripqp3_lp.csv"), DataFrame) 
ripqp_ma57_lp = CSV.read(string(res_path, "\\ripqp_ma571_lp.csv"), DataFrame)
ripqp_ma97_lp = CSV.read(string(res_path, "\\ripqp_ma971_lp.csv"), DataFrame)
ripqp_ma57_multi_lp = CSV.read(string(res_path, "\\ripqp_ma57_multi1_lp.csv"), DataFrame)
ripqp_ma57nosqd_lp = CSV.read(string(res_path, "\\ripqp_ma57nosqd1_lp.csv"), DataFrame)
ripqp_ma57_lp2 = CSV.read(string(res_path, "\\ripqp_ma572_lp.csv"), DataFrame) # regu comm dhab
ripqp_ma57nosqd_lp2 = CSV.read(string(res_path, "\\ripqp_ma57nosqd2_lp.csv"), DataFrame)
ripqp_qdldl_lp = CSV.read(string(res_path, "\\ripqp_qdldl1_lp.csv"), DataFrame)
ripqp_cholmod_lp = CSV.read(string(res_path, "\\ripqp_cholmod1_lp.csv"), DataFrame)
ripqp_multi1_lp = CSV.read(string(res_path, "\\ripqp_multi1_lp.csv"), DataFrame)
ripqp_nops1_lp = CSV.read(string(res_path, "\\ripqp_nops1_lp.csv"), DataFrame)
ripqp_cc1_lp = CSV.read(string(res_path, "\\ripqp_cc1_lp.csv"), DataFrame)

easy_pbs_lp = findall(ripqp1_lp.elapsed_time .≤ 0.5)
stats_lp = Dict(
                # :gurobi => gurobi1_lp,
                # :cplex => cplex1_lp,
                # :xpress => xpress1_lp,
                # :gurobi_nops1 => gurobi_nops1_lp,
                # :cplex_nops1 => cplex_nops1_lp,
                # :xpress_nops1 => xpress_nops1_lp,
                :ripqp => ripqp1_lp,
                # :ripqp2 => ripqp2_lp,
                # :ripqp3 => ripqp3_lp,
                # :ripqp_ldlfact => ripqp1_lp,
                # :ripqp_multi => ripqp_multi1_lp,
                # :ripqp_nops1 => ripqp_nops1_lp,
                #  :ripqp_ma57 => ripqp_ma57_lp,
                #  :ripqp_ma97 => ripqp_ma97_lp,
                #  :ripqp_ma57_multi => ripqp_ma57_multi_lp,
                # :ripqp_ma572 => ripqp_ma57_lp2,
                # :ripqp => filter(x -> x.id ∉ easy_pbs_lp, ripqp1_lp),
                # :ripqp_ma97 => filter(x -> x.id ∉ easy_pbs_lp, ripqp_ma97_lp),
                # :ripqp_ma57 => filter(x -> x.id ∉ easy_pbs_lp, ripqp_ma57_lp),
                # :ripqp_ma57nosqd => ripqp_ma57nosqd_lp2,
                # :ripqp_qdldl => ripqp_qdldl_lp,
                # :ripqp_cholmod => ripqp_cholmod_lp,
                :ripqp_cc => ripqp_cc1_lp,
                )

gurobi1_qp = CSV.read(string(res_path, "\\gurobi1_qp.csv"), DataFrame)
cplex1_qp = CSV.read(string(res_path, "\\cplex1_qp.csv"), DataFrame)
xpress1_qp = CSV.read(string(res_path, "\\xpress1_qp.csv"), DataFrame)
gurobi_nops1_qp = CSV.read(string(res_path, "\\gurobi_nops1_qp.csv"), DataFrame)
cplex_nops1_qp = CSV.read(string(res_path, "\\cplex_nops1_qp.csv"), DataFrame)
xpress_nops1_qp = CSV.read(string(res_path, "\\xpress_nops1_qp.csv"), DataFrame)
ripqp1_qp = CSV.read(string(res_path, "\\ripqp1_qp.csv"), DataFrame)
ripqp2_qp = CSV.read(string(res_path, "\\ripqp2_qp.csv"), DataFrame)
ripqp3_qp = CSV.read(string(res_path, "\\ripqp3_qp.csv"), DataFrame)
ripqp_ma57_qp = CSV.read(string(res_path, "\\ripqp_ma571_qp.csv"), DataFrame)
ripqp_ma97_qp = CSV.read(string(res_path, "\\ripqp_ma971_qp.csv"), DataFrame)
ripqp_ma57_multi_qp = CSV.read(string(res_path, "\\ripqp_ma57_multi1_qp.csv"), DataFrame)
ripqp_ma57nosqd_qp = CSV.read(string(res_path, "\\ripqp_ma57nosqd1_qp.csv"), DataFrame)
ripqp_ma57_qp2 = CSV.read(string(res_path, "\\ripqp_ma572_qp.csv"), DataFrame)
ripqp_ma57nosqd_qp2 = CSV.read(string(res_path, "\\ripqp_ma57nosqd2_qp.csv"), DataFrame)
ripqp_qdldl_qp = CSV.read(string(res_path, "\\ripqp_qdldl1_qp.csv"), DataFrame)
ripqp_cholmod_qp = CSV.read(string(res_path, "\\ripqp_cholmod1_qp.csv"), DataFrame)
ripqp_multi1_qp = CSV.read(string(res_path, "\\ripqp_multi1_qp.csv"), DataFrame)
ripqp_nops1_qp = CSV.read(string(res_path, "\\ripqp_nops1_qp.csv"), DataFrame)
ripqp_cc1_qp = CSV.read(string(res_path, "\\ripqp_cc1_qp.csv"), DataFrame)

easy_pbs_qp = findall(ripqp_ma57_qp.elapsed_time .≤ 1.0)
stats_qp = Dict(
                # :gurobi => gurobi1_qp,
                # :cplex => cplex1_qp,
                # :xpress => xpress1_qp,
                # :gurobi_nops1 => gurobi_nops1_qp,
                # :cplex_nops1 => cplex_nops1_qp,
                # :xpress_nops1 => xpress_nops1_qp,
                :ripqp => ripqp1_qp,
                # :ripqp2 => ripqp2_qp,
                # :ripqp_ma57 => ripqp_ma57_qp,
                # :ripqp_ma97 => ripqp_ma97_qp,
                # :ripqp_ma57_multi => ripqp_ma57_multi_qp,
                # :ripqp_ma572 => ripqp_ma57_qp2,
                # :ripqp_ma57nosqd => ripqp_ma57nosqd_qp,
                # :ripqp_qdldl => ripqp_qdldl_qp,
                # :ripqp => filter(x -> x.id ∉ easy_pbs_qp, ripqp_ma57_qp),
                # :ripqp_qdldl => filter(x -> x.id ∉ easy_pbs_qp, ripqp_qdldl_qp),
                # :ripqp_ma57_multi => filter(x -> x.id ∉ easy_pbs_qp, ripqp_ma57_multi_qp),
                # :ripqp_cholmod => ripqp_cholmod_qp,
                # :ripqp_multi => ripqp_multi1_qp,
                # :ripqp_nops1 => ripqp_nops1_qp,
                :ripqp_cc => ripqp_cc1_qp,
                )

function dfstat(df)
  output = zeros(length(df.status))
  for i=1:length(df.status)
    if df.primal_feas[i] === missing
      output[i] = Inf
    else 
      # output[i] = df.iter[i]
      # output[i] = df.relative_iter_cnt[i]
      output[i] = df.elapsed_time[i]
    end
    if df.status[i] ∉ ["first_order", "acceptable"]
      output[i] = Inf
    end
  end
  return output
end

cost = df -> df.elapsed_time + (df.status .!= :first_order) * Inf # + (df.elapsed_time .>= 10.) * Inf
pgfplotsx()
# perf = performance_profile(stats_lp, dfstat,legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# title!("Performance profile (Netlib problems)")
perf = performance_profile(stats_qp, dfstat,legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
title!("Performance profile (Maros and Meszaros problems)")
# display("image/svg+xml", perf)
# savefig(perf, raw"C:\Users\Geoffroy Leconte\Documents\doctorat\biblio\papiers\ripqp\paper\profiles\cc_mm_time.tikz")