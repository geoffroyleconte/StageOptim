using RipQPBenchmarks

save_path = "/home/gelecd/code/docGL/benchmarks/frontal22_results/ripqp_paper_bm"

ripqp_all_benchmarks(
  save_path;
  run_cplex = true,
  run_gurobi = true,
  run_xpress = true,
  run_ma57 = true,
  run_ma97 = true,
  plot_extension = ".pdf",
)
