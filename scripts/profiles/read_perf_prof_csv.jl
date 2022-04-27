using Plots, CSV, BenchmarkProfiles, DelimitedFiles, DataFrames, SolverBenchmark

# path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test2"
# path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test_qp2"
path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\benchmarks\frontal22_results\prof1_lp"
# solvers_list = Vector{String}(readdlm(string(path, "/test2_solvs.txt"))[:])
stats_solv = Dict{Symbol, DataFrame}()
for solver in readdir(path)[1:end-1][[9, 14]]
  # 3 , 4 , 9 , 13 (full reortho) , 14 , 17 , 21, 23 (full reortho), 24, 41 (full reortho) ,  
  ripqp_stats = CSV.read(string(path, "/", string(solver)), DataFrame)
  stats_solv[Symbol(solver[1:end-4])] = ripqp_stats
end
for i in 1:length(readdir(path))
  println((i, readdir(path)[i]))
end

function dfstat(df)
  output = zeros(length(df.pdd))
  for i=1:length(df.pdd)
    if df.pdd[i] === missing
      output[i] = Inf
    else 
      output[i] = df.absolute_iter_cnt[i]
      # output[i] = df.elapsed_time[i]
    end
    if df.status[i] != "acceptable"
      output[i] = Inf
    end
  end
  return output
end

perf = performance_profile(stats_solv, dfstat, legend=:topright)
# ylims!((0.0, 0.9))
# plot(perf, )
title!("time perf profile (Netlib problems)")
# title!("iterations perf profile (Netlib problems)")
# title!("time perf profile (Maros and Meszaros problems)")
# title!("pdd perf profile (Netlib problems)")

display("image/svg+xml", perf)
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\biblio\papiers\proposition_recherche\prop_rech_GL\images\perf_ksysts"
# savefig(string(save_path, "/normal2.pdf"))

# select(stats_solv[:K2_LDL], [2, 3, 4, 6, 7, 8, 10, 11, 40, 41])
# df[!,:status] = convert.(String,df[!,:status])
# df[!,:name] = convert.(String,df[!,:name])
# table = pretty_table(df, backend = Val(:latex))
# latexify(df)