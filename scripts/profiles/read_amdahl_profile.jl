using Plots, JLD, BenchmarkProfiles

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\data_profiles"
data_solv = load(string(path, "/test2.jld"))["data"]
solver_list = readdlm(string(path, "/test2_solvs.txt"))
n_iter, n_pb, n_solvers = size(data_solv)
N = ones(n_pb)

perf = data_profile(PlotsBackend(), data_solv, N, solvers_list, legend=:topright,
                    τ= 1.0e-3)
# plot(perf, )
title!("data profile (Netlib problems)")

display("image/svg+xml", perf)
