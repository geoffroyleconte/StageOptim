using Plots, JLD, BenchmarkProfiles, DelimitedFiles

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\data_profiles"
n_iter, n_pb, n_solvers = size(data_solv)
N = ones(n_pb)
data_solv = load(string(path, "/test2.jld"))["data"]
solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])

perf = data_profile(PlotsBackend(), data_solv, N, solvers_list, legend=:topleft, τ= 1.0e-3)
# plot(perf, )
title!("data profile (Netlib problems)")

display("image/svg+xml", perf)
