using Plots, JLD, BenchmarkProfiles, DelimitedFiles

path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\data_profiles"
data_solv = load(string(path, "/test2.jld"))["data"]
n_iter, n_pb, n_solvers = size(data_solv)
n_pb = findfirst(isequal(0.0), data_solv[end, :, 1])-1
N = ones(n_pb)
solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
println(solvers_list)

data_solv_r = data_solv[:, 1:n_pb, [4, end-3, end-1]]
solvers_list_r = solvers_list[[4, end-3, end-1]]
perf = data_profile(PlotsBackend(), data_solv_r, N, solvers_list_r, 
                    legend=:topleft, τ= 1.0e-3)
ylims!((0.0, 0.9))
# plot(perf, )
title!("data profile (Netlib problems)")

display("image/svg+xml", perf)
