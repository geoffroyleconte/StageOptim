# res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\amdahl_benchmarks\\results"
res_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\benchmarks\ripqp_paper"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\biblio\papiers\ripqp\paper\profiles"
using Plots
using DataFrames, SolverBenchmark, SolverTools
using JLD2
using CSV
# using FileIO
open_file(fname; res_path = res_path) = CSV.read(string(res_path, "\\", fname, ".csv"), DataFrame)

gurobi1_lp = open_file("gurobi1_lp")
cplex1_lp = open_file("cplex1_lp")
xpress1_lp = open_file("xpress1_lp")
ripqp1_lp = open_file("ripqp1_lp") # compare commercial solvs + cc1
ripqp_cc1_lp = open_file("ripqp_cc1_lp")
ripqp_ma57_lp = open_file("ripqp_ma57_lp")
ripqp_ma97_lp = open_file("ripqp_ma971_lp")
ripqp_qdldl_lp = open_file("ripqp_qdldl1_lp")
ripqp_cholmod_lp = open_file("ripqp_cholmod1_lp")
ripqp_multi1_lp = open_file("ripqp_multi1_lp")
ripqp_ma57_multi1_lp = open_file("ripqp_ma57_multi1_lp")
ripqp_ma57_multi2_lp = open_file("ripqp_ma57_multi2_lp")
ripqp_ldlprecond1_lp = open_file("ripqp_ldlprecond1_lp") # regu1 1.0e-8, stop crit 64, no equi
ripqp_ldlprecond2_lp = open_file("ripqp_ldlprecond2_lp") # regu1 1.0e-8 equi
ripqp_lldlprecond1_lp = open_file("ripqp_lldlprecond_lp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch

easy_pbs_lp = findall(ripqp1_lp.elapsed_time .≤ 10.0)

gurobi1_qp = open_file("gurobi1_qp")
cplex1_qp = open_file("cplex1_qp")
xpress1_qp = open_file("xpress1_qp")
ripqp1_qp = open_file("ripqp1_qp") # compare commercial solvs + cc1
ripqp_cc1_qp = open_file("ripqp_cc1_qp")
ripqp_ma57_qp = open_file("ripqp_ma57_qp")
ripqp_ma97_qp = open_file("ripqp_ma971_qp")
ripqp_qdldl_qp = open_file("ripqp_qdldl1_qp")
ripqp_cholmod_qp = open_file("ripqp_cholmod1_qp")
ripqp_multi1_qp = open_file("ripqp_multi1_qp")
ripqp_ma57_multi1_qp = open_file("ripqp_ma57_multi1_qp")
ripqp_ma57_multi2_qp = open_file("ripqp_ma57_multi2_qp")
ripqp_ldlprecond1_qp = open_file("ripqp_ldlprecond1_qp") # regu1 1.0e-8, stop crit 64, no equi
ripqp_ldlprecond2_qp = open_file("ripqp_ldlprecond2_qp") # regu1 1.0e-8 equi
ripqp_lldlprecond1_qp = open_file("ripqp_lldlprecond_qp") # regu1 1.0e-8, stop crit 64, no equi, new regu_try_catch

function dfstat_time(df)
  output = zeros(length(df.status))
  for i=1:length(df.status)
    if df.primal_feas[i] === missing || df.objective[i] == Inf
      output[i] = Inf
    else 
      # output[i] = df.iter[i]
      # output[i] = df.relative_iter_cnt[i]
      # output[i] = df.iters_sp[i]
      output[i] = df.elapsed_time[i]
      # output[i] = 4 * df.iters_sp2[i] + df.iters_sp[i]
    end
    if df.status[i] ∉ ["first_order", "acceptable"]
      output[i] = Inf
    end
  end
  return output
end

function dfstat_energy(df)
  output = zeros(length(df.status))
  for i=1:length(df.status)
    if df.primal_feas[i] === missing || df.objective[i] == Inf
      output[i] = Inf
    else 
      output[i] = df.relative_iter_cnt[i]
      # output[i] = 4 * df.iters_sp2[i] + df.iters_sp[i]
    end
    if df.status[i] ∉ ["first_order", "acceptable"]
      output[i] = Inf
    end
  end
  return output
end

function dfstat_energy2(df)
  output = zeros(length(df.status))
  for i=1:length(df.status)
    if df.primal_feas[i] === missing || df.objective[i] == Inf
      output[i] = Inf
    else 
      output[i] = 4 * df.iters_sp2[i] + df.iters_sp[i]
    end
    if df.status[i] ∉ ["first_order", "acceptable"]
      output[i] = Inf
    end
  end
  return output
end

cost = df -> df.elapsed_time + (df.status .!= :first_order) * Inf # + (df.elapsed_time .>= 10.) * Inf

################################ compare solvers ################################
# stats_lp = Dict(:gurobi => gurobi1_lp, :cplex => cplex1_lp, :xpress => xpress1_lp, :ripqp => ripqp1_lp)
# stats_qp = Dict(:gurobi => gurobi1_qp, :cplex => cplex1_qp, :xpress => xpress1_qp, :ripqp => ripqp1_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\solvers_net_time.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\solvers_mm_time.tikz"))

# ################################ compare fact alg ################################
# stats_lp = Dict(:ripqp_ldl => ripqp1_lp, :ripqp_ma57 => ripqp_ma57_lp, :ripqp_ma97 => ripqp_ma97_lp, :ripqp_qdldl => ripqp_qdldl_lp, :ripqp_cholmod => ripqp_cholmod_lp)
# stats_qp = Dict(:ripqp_ldl => ripqp1_qp, :ripqp_ma57 => ripqp_ma57_qp, :ripqp_ma97 => ripqp_ma97_qp, :ripqp_qdldl => ripqp_qdldl_qp, :ripqp_cholmod => ripqp_cholmod_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\fact_net_time.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\fact_mm_time.tikz"))

# ################################ centrality corrections ################################
# stats_lp = Dict(:ripqp => ripqp1_lp, :ripqp_cc => ripqp_cc1_lp)
# stats_qp = Dict(:ripqp => ripqp1_qp, :ripqp_cc => ripqp_cc1_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\cc_net_time.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\cc_mm_time.tikz"))

# ################################ multi precision ################################
# stats_lp = Dict(:ripqp_mono => ripqp1_lp, :ripqp_multi => ripqp_multi1_lp)
# stats_qp = Dict(:ripqp_mono => ripqp1_qp, :ripqp_multi => ripqp_multi1_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_energy, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\multi_net_riter.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_energy, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\multi_mm_riter.tikz"))

# ################################ multi precision ma57 ################################
# stats_lp = Dict(:ripqp_ma57_mono => ripqp_ma57_lp, :ripqp_ma57_multi1 => ripqp_ma57_multi1_lp, :ripqp_ma57_multi2 => ripqp_ma57_multi2_lp)
# stats_qp = Dict(:ripqp_ma57_mono => ripqp_ma57_qp, :ripqp_ma57_multi1 => ripqp_ma57_multi1_qp, :ripqp_ma57_multi2 => ripqp_ma57_multi2_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\ma57_net_time.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\ma57_mm_time.tikz"))

# ################################ multi precision ldl precond ################################
# stats_lp = Dict(:ripqp_multi => ripqp_multi1_lp, :ripqp_multifact1 => ripqp_ldlprecond1_lp, :ripqp_multifact2 => ripqp_ldlprecond2_lp)
# stats_qp = Dict(:ripqp_multi => ripqp_multi1_qp, :ripqp_multifact1 => ripqp_ldlprecond1_qp, :ripqp_multifact2 => ripqp_ldlprecond2_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_energy2, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\multifact_net_riter.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_energy2, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\multifact_mm_riter.tikz"))

# ################################ multi precision lldl precond ################################
# stats_lp = Dict(:ripqp_multifact1 => ripqp_ldlprecond1_lp, :ripqp_multifact_lldl => ripqp_lldlprecond1_lp)
# stats_qp = Dict(:ripqp_multifact1 => ripqp_ldlprecond1_qp, :ripqp_multifact_lldl => ripqp_lldlprecond1_qp)
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_energy2, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\lldl_net_riter.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_energy2, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\lldl_mm_riter.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_lp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\lldl_net_time.tikz"))
# pgfplotsx()
# perf = performance_profile(stats_qp, dfstat_time, legend=:bottomright, b = SolverBenchmark.BenchmarkProfiles.PGFPlotsXBackend())
# # savefig(perf, string(save_path, "\\lldl_mm_time.tikz"))


save_tbl_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\biblio\papiers\ripqp\paper\tables"
using PrettyTables

safe_latex_customstring(s::AbstractString) = "\\(" * replace(s, "_" => "\\_") * "\\)"
function safe_latex_customstring(col::Integer)
  # by this point, the table value should already have been converted to a string
  return (s, i, j) -> begin
    j != col && return s
    return safe_latex_customstring(s)
  end
end

# SolverBenchmark.default_formatters[InlineString] = "%15s"
hdr_override = Dict(
  :name => "Name",
  :nvar => "n",
  :ncon => "m",
  :elapsed_time => "time (s)",
  :iter => "iter tot",
  :iters_sp => "iter32",
  :iters_sp2 => "iter64",
  :primal_feas => "pfeas",
  :dual_feas => "dfeas",
)

fmt_override = Dict(
  :status => "15s",
  :objective => "%8.1e",
  :elapsed_time => "%8.1e",
  :primal_feas => "%8.1e",
  :dual_feas => "%8.1e",
  :pdd => "%8.1e",
)
pretty_stats(stdout,
  ripqp1_lp[!, [:name, :nvar, :ncon, :status, :objective, :pdd, :primal_feas, :dual_feas, :elapsed_time, :iter, :iters_sp, :iters_sp2]],
  hdr_override = hdr_override,
  col_formatters = fmt_override,)

open(string(save_tbl_path, "\\ripqp_multi1_qp.tex"), "w") do io
  pretty_latex_stats(io,
    ripqp_multi1_qp[!, [:name, :nvar, :ncon, :status, :objective, :pdd, :primal_feas, :dual_feas, :elapsed_time, :iter, :iters_sp, :iters_sp2]],
    hdr_override = hdr_override,
    col_formatters = fmt_override,
    formatters = (ft_printf(["15s", "%7.1e", "%7.1e","%7.1e","%7.1e", "%7.1e"], [4, 5, 6,7,8, 9]),
    SolverBenchmark.safe_latex_AbstractFloat_col(5),
    SolverBenchmark.safe_latex_AbstractFloat_col(6),
    SolverBenchmark.safe_latex_AbstractFloat_col(7),
    SolverBenchmark.safe_latex_AbstractFloat_col(8),
    SolverBenchmark.safe_latex_AbstractFloat_col(9),
    safe_latex_customstring(4)))
end