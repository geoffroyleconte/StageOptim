using Plots, JLD2, BenchmarkProfiles, DelimitedFiles, DataFrames, SolverBenchmark

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test2"
# path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test_qp1"
solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
stats_solv = Dict{Symbol, DataFrame}()
for solver in solvers_list[[2, 4, 6, end-3, end-2]]
  file = jldopen(string(path, "/", string(solver), ".jld2"), "r")
  ripqp_stats = file["stats"];
  close(file)
  stats_solv[Symbol(solver)] = ripqp_stats
end
println(solvers_list)

function dfpdd(df)
  output = zeros(length(df.pdd))
  for i=1:length(df.pdd)
    if df.pdd[i] === missing
      output[i] = Inf
    else 
      output[i] = df.absolute_iter_cnt[i]
    end
    if df.status[i] != :acceptable
      output[i] = Inf
    end
  end
  return output
end

perf = performance_profile(stats_solv, dfpdd, legend=:bottomright)
# ylims!((0.0, 0.9))
# plot(perf, )
title!("iter performance profile (Netlib problems)")
# title!("pdd perf profile (Netlib problems)")

display("image/svg+xml", perf)
