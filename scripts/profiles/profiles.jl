# res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\frontal22_results\\bm1"
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

# rip_mono = open_file(string(res_path, "\\G-2021-03_lp_mono.jld2"));
# rip_mono1 = open_file(string(res_path, "\\ripqp_mono_sc1_lp.jld2"))#, "r");
# rip_mono2 = open_file(string(res_path, "\\ripqp_mono_sc2_lp.jld2"));
rip_mono1 = CSV.read(string(res_path, "\\ripqp_mono_classic_lp.csv"), DataFrame)
rip_mono2 = CSV.read(string(res_path, "\\ripqp_mono_sc2_lp.csv"), DataFrame)
rip_mono3 = CSV.read(string(res_path, "\\ripqp_mono_wtimerout_lp.csv"), DataFrame)
# rip_multiK2 = open_file(string(res_path, "\\ripqp_multi_K2_lp.jld2"));
# rip_mono_c = open_file(string(res_path, "\\ripqp_ccorr_1_lp.jld2")); 
# rip_multi_z = open_file(string(res_path, "\\ripqp_multi_z_lp.jld2"))
# rip_multi_r = open_file(string(res_path, "\\ripqp_multi_r_lp.jld2"))
# rip_multi1 = open_file(string(res_path, "\\ripqp_multi_1_lp.jld2"));
# rip_minres1 = open_file(string(res_path, "\\ripqp_minres_1_lp.jld2"));
# rip_minres2 = open_file(string(res_path, "\\ripqp_minres_2_lp.jld2"));
# rip_monoIPFK22 = open_file(string(res_path, "\\ripqp_mono_IPFK2_2_lp.jld2")); # r, γ =  T(0.95), T(0.1)
# rip_monoIPFK23 = open_file(string(res_path, "\\ripqp_mono_IPFK2_3_lp.jld2")); # r, γ =  T(0.999), T(0.05), fix bug

stats_lp = Dict(
                :ripqp_classic => rip_mono1,
                # :ripqp_newscale => rip_mono2,
                :ripqp_to => rip_mono3,
                # :rip_monoIPFK22 => rip_monoIPFK22,
                # :rip_monoIPFK23 => rip_monoIPFK23,
                )

cost = df -> df.elapsed_time + (df.status .== :acceptable) * Inf# + (df.elapsed_time .>= 10.) * Inf
perf = performance_profile(stats_lp, cost,legend=:bottomright)
title!("Performance profile (Maros Meszaros problems)")
# savefig(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\profiles\minres_net_iter_wrip.pdf")

################################ QP #####################
# rip_mono = open_file(string(res_path, "\\G-2021-03_qp_mono.jld2"));
# rip_mono1 = open_file(string(res_path, "\\ripqp_mono_1_qp.jld2"));
# rip_multiK2 = open_file(string(res_path, "\\ripqp_multi_K2_qp.jld2"));
# rip_mono_c = open_file(string(res_path, "\\ripqp_ccorr_1_qp.jld2")); 
# rip_multi_z = open_file(string(res_path, "\\ripqp_multi_z_qp.jld2"))
# rip_multi_r = open_file(string(res_path, "\\ripqp_multi_r_qp.jld2"))
# rip_multi1 = open_file(string(res_path, "\\ripqp_multi_1_qp.jld2"));
# rip_multi_d1 = open_file(string(res_path, "\\ripqp_multi_dynamic_1_qp.jld2"));

# stats_qp = Dict(
#                 :ripqp_mono1 => rip_mono1,
#                 :ripqp_multiK2 => rip_multiK2,
#                 :ripqp_multi_r => rip_multi_r,
#                 :ripqp_multi_z => rip_multi_z,
#                 )

# perf = performance_profile(stats_qp, df->df.iter)
# plot!(perf, legend=:bottomright)
# title!("Performance profile (Maros and Meszaros problems)")
# # savefig(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\profiles\ccorr_mmwmono.pdf")
