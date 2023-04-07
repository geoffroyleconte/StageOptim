using RipQPBenchmarks

save_path = "/home/gelecd/code/docGL/benchmarks/frontal22_results/ripqp_paper_bm"
run_benchmarks_solvers(
  save_path;
  run_cplex = true,
  run_gurobi = true,
  run_xpress = false,
  run_ma57 = true,
  run_ma97 = true,
)

run_benchmarks_quad(save_path)
