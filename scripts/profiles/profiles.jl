res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
using DataFrames, SolverBenchmark, SolverTools, JLD2, Plots

function open_file(path)
    file = jldopen(path, "r")
    df_out = file["stats"]
    close(file)
    return df_out
end

rip_mono = open_file(string(res_path, "\\G-2021-03_lp_mono.jld2"));
rip_mono1 = open_file(string(res_path, "\\ripqp_mono_1_lp.jld2"));
rip_multiK2 = open_file(string(res_path, "\\ripqp_multi_K2_lp.jld2"));
rip_mono_c = open_file(string(res_path, "\\ripqp_ccorr_1_lp.jld2")); 
rip_multi_z = open_file(string(res_path, "\\ripqp_multi_z_lp.jld2"))
rip_multi_r = open_file(string(res_path, "\\ripqp_multi_r_lp.jld2"))
rip_multi1 = open_file(string(res_path, "\\ripqp_multi_1_lp.jld2"));
rip_minres1 = open_file(string(res_path, "\\ripqp_minres_1_lp.jld2"));
rip_minres2 = open_file(string(res_path, "\\ripqp_minres_2_lp.jld2"));

stats_lp = Dict(
                :ripqp_mono1 => rip_mono1,
                :ripqp_minres1 => rip_minres1,
                :ripqp_minres2 => rip_minres2,
                )

cost = df -> df.iter + (df.status .!= :acceptable) * Inf + (df.elapsed_time .>= 10.) * Inf
perf = performance_profile(stats_lp, cost,legend=:topright)
title!("Performance profile (Netlib problems)")
# savefig(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\profiles\minres_net_iter_wrip.pdf")

################################ QP #####################
rip_mono = open_file(string(res_path, "\\G-2021-03_qp_mono.jld2"));
rip_mono1 = open_file(string(res_path, "\\ripqp_mono_1_qp.jld2"));
rip_multiK2 = open_file(string(res_path, "\\ripqp_multi_K2_qp.jld2"));
rip_mono_c = open_file(string(res_path, "\\ripqp_ccorr_1_qp.jld2")); 
rip_multi_z = open_file(string(res_path, "\\ripqp_multi_z_qp.jld2"))
rip_multi_r = open_file(string(res_path, "\\ripqp_multi_r_qp.jld2"))
rip_multi1 = open_file(string(res_path, "\\ripqp_multi_1_qp.jld2"));
rip_multi_d1 = open_file(string(res_path, "\\ripqp_multi_dynamic_1_qp.jld2"));

stats_qp = Dict(
                :ripqp_mono1 => rip_mono1,
                :ripqp_multiK2 => rip_multiK2,
                :ripqp_multi_r => rip_multi_r,
                :ripqp_multi_z => rip_multi_z,
                )

perf = performance_profile(stats_qp, df->df.iter)
plot!(perf, legend=:bottomright)
title!("Performance profile (Maros and Meszaros problems)")
# savefig(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\profiles\ccorr_mmwmono.pdf")
