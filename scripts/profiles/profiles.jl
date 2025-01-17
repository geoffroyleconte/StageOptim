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
rip_mono = CSV.read(string(res_path, "\\ripqp_mono2_lp.csv"), DataFrame)
rip_multi = CSV.read(string(res_path, "\\ripqp_multi2_lp.csv"), DataFrame)
rip_multiref2 = CSV.read(string(res_path, "\\ripqp_multiref2_qp.csv"), DataFrame) # multi old v
rip_multizoom2 = CSV.read(string(res_path, "\\ripqp_multizoom2_qp.csv"), DataFrame) # multi old v
rip_multiref4 = CSV.read(string(res_path, "\\ripqp_multiref4_qp.csv"), DataFrame) # old v start pts
rip_multizoom4 = CSV.read(string(res_path, "\\ripqp_multizoom4_qp.csv"), DataFrame)
rip_multiref5 = CSV.read(string(res_path, "\\ripqp_multiref5_qp.csv"), DataFrame) # old v start pts + s_l
rip_multizoom5 = CSV.read(string(res_path, "\\ripqp_multizoom5_qp.csv"), DataFrame) # 
rip_multiref6 = CSV.read(string(res_path, "\\ripqp_multiref6_qp.csv"), DataFrame) # new w start pt
rip_multizoom6 = CSV.read(string(res_path, "\\ripqp_multizoom6_lp.csv"), DataFrame) # 
rip_multiref7 = CSV.read(string(res_path, "\\ripqp_multiref7_qp.csv"), DataFrame) # 
rip_multizoom7 = CSV.read(string(res_path, "\\ripqp_multizoom7_lp.csv"), DataFrame) # 
# ripqp_eq = CSV.read(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\benchmarks\frontal22_results\bm1\ripqp_equi_minres_qlp_lp.csv", DataFrame)
# rip_multiK2 = open_file(string(res_path, "\\ripqp_multi_K2_lp.jld2"));
# rip_mono_c = open_file(string(res_path, "\\ripqp_ccorr_1_lp.jld2")); 
# rip_multi1 = open_file(string(res_path, "\\ripqp_multi_1_lp.jld2"));
# rip_minres1 = open_file(string(res_path, "\\ripqp_minres_1_lp.jld2"));
# rip_minres2 = open_file(string(res_path, "\\ripqp_minres_2_lp.jld2"));
# rip_monoIPFK22 = open_file(string(res_path, "\\ripqp_mono_IPFK2_2_lp.jld2")); # r, γ =  T(0.95), T(0.1)
# rip_monoIPFK23 = open_file(string(res_path, "\\ripqp_mono_IPFK2_3_lp.jld2")); # r, γ =  T(0.999), T(0.05), fix bug

stats_lp = Dict(
                # :ripqp_classic => rip_mono1,
                # :ripqp_newscale => rip_mono2,
                # :mono => rip_mono,
                :multi => rip_multi,
                # :multiref => rip_multiref,
                # :multizoom => rip_multizoom,
                # :multiref2 => rip_multiref2,
                # :multizoom2 => rip_multizoom2,
                # :multizoom4 => rip_multizoom4,
                # :multizoom5 => rip_multizoom5,
                :multizoom => rip_multizoom6,
                :multizoom7 => rip_multizoom7,
                # :multiref6 => rip_multiref6,
                # :multiref5 => rip_multiref5,
                # :multiref4 => rip_multiref4,
                # :multizoom6 => rip_multizoom6,
                # :ripqp2 => rip_mono2,
                # :K2Equilibration_minresqlp => ripqp_eq,
                # :rip_monoIPFK22 => rip_monoIPFK22,
                # :rip_monoIPFK23 => rip_monoIPFK23,
                )

function dfstat(df)
  output = zeros(length(df.pdd))
  for i=1:length(df.pdd)
    if df.pdd[i] === missing
      output[i] = Inf
    else 
    #   output[i] = df.absolute_iter_cnt[i]
      output[i] = df.relative_iter_cnt[i]
    end
    if df.status[i] != "first_order"
      output[i] = Inf
    end
  end
  return output
end

# cost = df -> df.elapsed_time + (df.status .!= :first_order) * Inf # + (df.elapsed_time .>= 10.) * Inf
perf = performance_profile(stats_lp, dfstat,legend=:bottomright)
title!("Performance profile (Netlib problems)")
# title!("Performance profile (Maros and Meszaros problems)")
display("image/svg+xml", perf)
# savefig(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\profiles\multizoom_iter_lp.pdf")

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
