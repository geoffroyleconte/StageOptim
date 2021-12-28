using Plots, JLD2, BenchmarkProfiles, DelimitedFiles, DataFrames, SolverBenchmark

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test1"
# path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test_qp1"
solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
stats_solv = Dict{Symbol, DataFrame}()
for solver in solvers_list[[1, 2, end-2, end-1, end]]
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
      output[i] = df.elapsed_time[i]
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
title!("time performance profile (Netlib problems)")
# title!("time perf profile (Maros and Meszaros problems)")

display("image/svg+xml", perf)
