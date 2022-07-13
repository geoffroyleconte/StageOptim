using Plots, CSV, BenchmarkProfiles, DelimitedFiles, DataFrames, SolverBenchmark

# path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test2"
# path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test_qp2"
path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\benchmarks\frontal22_results\prof2_lp"
# solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
stats_solv = Dict{Symbol, DataFrame}()
for solver in readdir(path)[1:end-1][[19,30, 40, 28, 13, 14]]
  ripqp_stats = CSV.read(string(path, "/", string(solver)), DataFrame)
  stats_solv[Symbol(solver[1:end-4])] = ripqp_stats
end
# stats_solv[:K2_minresqlp_equi_full] = CSV.read(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\benchmarks\frontal22_results\bm1\ripqp_equi_minres_qlp_lp.csv", DataFrame)
for i in 1:length(readdir(path))
  println((i, readdir(path)[i]))
end
# println((length(readdir(path))+1, "ripqp_equi_minres_qlp_lp.csv"))

function dfstat(df)
  output = zeros(length(df.pdd))
  for i=1:length(df.pdd)
    if df.pdd[i] === missing
      output[i] = Inf
    else 
      output[i] = df.absolute_iter_cnt[i]
      # output[i] = df.elapsed_time[i]
    end
    if df.status[i] != "first_order"
      output[i] = Inf
    end
  end
  return output
end

perf = performance_profile(stats_solv, dfstat, legend=:topright)
# ylims!((0.0, 0.9))
# plot(perf, )
title!("iter perf profile (Netlib problems)")
# title!("iterations perf profile (Netlib problems)")
# title!("time perf profile (Maros and Meszaros problems)")
# title!("pdd perf profile (Netlib problems)")

display("image/svg+xml", perf)
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\biblio\resume\rapports\2022-06-27"
# savefig(string(save_path, "/k1.pdf"))

# select(stats_solv[:K2_LDL], [2, 3, 4, 6, 7, 8, 10, 11, 40, 41])
# df[!,:status] = convert.(String,df[!,:status])
# df[!,:name] = convert.(String,df[!,:name])
# table = pretty_table(df, backend = Val(:latex))
# latexify(df)