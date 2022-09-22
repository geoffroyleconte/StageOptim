res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\ripqp_paper"
using Plots
using DataFrames, SolverBenchmark, SolverTools
using JLD2
using CSV
# using FileIO
open_file(fname; res_path = res_path) = CSV.read(string(res_path, "\\", fname, ".csv"), DataFrame)

ripqp_multiks = [
  open_file("ripqp_multik1_quad"), # 3 sp init 
  open_file("ripqp_multik2_quad"), # 
  open_file("ripqp_multik3_quad"), # 1 100iter gmres ir no eq
  open_file("ripqp_multik4_quad"), # 
  open_file("ripqp_multik5_quad"), # 1 100 iter gmres no eq
  open_file("ripqp_multik6_quad"), # 
  open_file("ripqp_multik7_quad"), #
] 
ripqp_multi1 = open_file("ripqp_multi1_quad")

nks = length(ripqp_multiks)
header = [
  "time",
  "iters",
  "iters 64",
  "iters 128",
  "objective",
  "pdd",
  "primal feas",
  "dual feas",
]
nh = length(header)
data = Matrix{Any}(undef, nks+1, nh)
solver_names = [string("k", i) for i in 1:nks]
push!(solver_names, "multi")
pb_index = 3
for i in 1:nks
  data[i, :] .= [
    ripqp_multiks[i].elapsed_time[pb_index],
    ripqp_multiks[i].iter[pb_index],
    ripqp_multiks[i].iters_sp[pb_index],
    ripqp_multiks[i].iters_sp2[pb_index],
    ripqp_multiks[i].objective[pb_index],
    ripqp_multiks[i].pdd[pb_index],
    ripqp_multiks[i].primal_feas[pb_index],
    ripqp_multiks[i].dual_feas[pb_index],
  ]
end
data[end, :] .= [
  ripqp_multi1.elapsed_time[pb_index],
  ripqp_multi1.iter[pb_index],
  ripqp_multi1.iters_sp[pb_index],
  ripqp_multi1.iters_sp2[pb_index],
  ripqp_multi1.objective[pb_index],
  ripqp_multi1.pdd[pb_index],
  ripqp_multi1.primal_feas[pb_index],
  ripqp_multi1.dual_feas[pb_index],
]

using PrettyTables


nh = length(header)
pbs = ["TMA_ME", "GlcAlift", "GlcAerWT"]

pretty_table(data; 
  header = header,
  row_names= solver_names,
  title = pbs[pb_index],
  # backend = Val(:latex),
  )