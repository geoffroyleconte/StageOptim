using Plots, JLD2, BenchmarkProfiles, DelimitedFiles, DataFrames, SolverBenchmark

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test1"
solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
stats_solv = Dict{Symbol, DataFrame}()
for solver in solvers_list
  file = jldopen(string(path, "/", string(solver), ".jld2"), "r")
  ripqp_stats = file["stats"];
  close(file)
  stats_solv[Symbol(solver)] = ripqp_stats
end
println(solvers_list)

perf = performance_profile(stats_solv, df->df.elapsed_time)
# ylims!((0.0, 0.9))
# plot(perf, )
title!("data profile (Netlib problems)")

display("image/svg+xml", perf)
