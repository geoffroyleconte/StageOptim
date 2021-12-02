using Plots, JLD, BenchmarkProfiles, DelimitedFiles

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\data_profiles"
data_solv = load(string(path, "/test2.jld"))["data"]
n_iter, n_pb, n_solvers = size(data_solv)
N = ones(n_pb)
solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
println(solvers_list)

data_solv_r = data_solv[:, :, 1:end]
solvers_list_r = solvers_list[1:end]
perf = data_profile(PlotsBackend(), data_solv_r, N, solvers_list_r, 
                    legend=:bottomright, τ= 1.0e-8)
ylims!((0.7, 1.0))
# plot(perf, )
title!("data profile (Netlib problems)")

display("image/svg+xml", perf)
