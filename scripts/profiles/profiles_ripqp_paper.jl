# res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\ripqp_paper"
using DataFrames, SolverBenchmark, SolverTools, Plots
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

stats_lp = Dict(
                :gurobi1 => gurobi1_lp,
                :cplex1 => cplex1_lp,
                :xpress1 => xpress1_lp,
                )

gurobi1_qp = CSV.read(string(res_path, "\\gurobi1_qp.csv"), DataFrame)
cplex1_qp = CSV.read(string(res_path, "\\cplex1_qp.csv"), DataFrame)

stats_qp = Dict(
                :gurobi1 => gurobi1_qp,
                :cplex1 => cplex1_qp,
                )

function dfstat(df)
  output = zeros(length(df.status))
  for i=1:length(df.status)
    if df.primal_feas[i] === missing
      output[i] = Inf
    else 
    #   output[i] = df.absolute_iter_cnt[i]
      output[i] = df.elapsed_time[i]
    end
    if df.status[i] ∉ ["first_order", "acceptable"]
      output[i] = Inf
    end
  end
  return output
end

# cost = df -> df.elapsed_time + (df.status .!= :first_order) * Inf # + (df.elapsed_time .>= 10.) * Inf
# perf = performance_profile(stats_lp, dfstat,legend=:bottomright)
# title!("Performance profile (Netlib problems)")
perf = performance_profile(stats_qp, dfstat,legend=:bottomright)
title!("Performance profile (Maros and Meszaros problems)")
display("image/svg+xml", perf)
# savefig(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\profiles\multizoom_iter_lp.pdf")
