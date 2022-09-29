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
ripqp_mono1 = open_file("ripqp_mono1_quad")
ripqp_multi_nops3 = open_file("ripqp_multik3_nops_quad") # en fait c'est le 4
ripqp_multi_nops5 = open_file("ripqp_multik5_nops_quad")

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
data = Matrix{Any}(undef, nks+2, nh)
solver_names = [string("k", i) for i in 1:nks]
push!(solver_names, "multi")
push!(solver_names, "mono")
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
data[end-1, :] .= [
  ripqp_multi1.elapsed_time[pb_index],
  ripqp_multi1.iter[pb_index],
  ripqp_multi1.iters_sp[pb_index],
  ripqp_multi1.iters_sp2[pb_index],
  ripqp_multi1.objective[pb_index],
  ripqp_multi1.pdd[pb_index],
  ripqp_multi1.primal_feas[pb_index],
  ripqp_multi1.dual_feas[pb_index],
]

data[end, :] .= [
  ripqp_mono1.elapsed_time[pb_index],
  ripqp_mono1.iter[pb_index],
  0,
  ripqp_mono1.iters_sp[pb_index],
  ripqp_mono1.objective[pb_index],
  ripqp_mono1.pdd[pb_index],
  ripqp_mono1.primal_feas[pb_index],
  ripqp_mono1.dual_feas[pb_index],
]

using PrettyTables


nh = length(header)
pbs = ["TMA ME", "GlcAlift", "GlcAerWT"]

# pretty_table(data; 
#   header = header,
#   row_names= solver_names,
#   title = pbs[pb_index],
#   # backend = Val(:latex),
#   formatters = ft_printf(["%7.1e", "%d", "%d", "%d", "%7.1e","%7.1e","%7.1e","%7.1e"], 1:8),
#   )
header2 = [
  "solver",
  "time",
  "iter64",
  "iter128",
  "obj",
  "pdd",
  "pfeas",
  "dfeas",
]
data2 = Matrix{Any}(undef, 12, nh)
row_names2 = []
for pb_index in 1:3
  push!(row_names2, pbs[pb_index])
  push!(row_names2, pbs[pb_index])
  push!(row_names2, pbs[pb_index])
  push!(row_names2, pbs[pb_index])
  data2[4 * (pb_index-1) + 1, :] .= [
    "multik",
    ripqp_multiks[3].elapsed_time[pb_index],
    ripqp_multiks[3].iters_sp[pb_index],
    ripqp_multiks[3].iters_sp2[pb_index],
    ripqp_multiks[3].objective[pb_index],
    ripqp_multiks[3].pdd[pb_index],
    ripqp_multiks[3].primal_feas[pb_index],
    ripqp_multiks[3].dual_feas[pb_index],
  ]
  data2[4 * (pb_index-1) + 2, :] .= [
    "multi",
    ripqp_multi1.elapsed_time[pb_index],
    ripqp_multi1.iters_sp[pb_index],
    ripqp_multi1.iters_sp2[pb_index],
    ripqp_multi1.objective[pb_index],
    ripqp_multi1.pdd[pb_index],
    ripqp_multi1.primal_feas[pb_index],
    ripqp_multi1.dual_feas[pb_index],
  ]

  data2[4 * (pb_index-1) + 3, :] .= [
    "mono",
    ripqp_mono1.elapsed_time[pb_index],
    0,
    ripqp_mono1.iters_sp[pb_index],
    ripqp_mono1.objective[pb_index],
    ripqp_mono1.pdd[pb_index],
    ripqp_mono1.primal_feas[pb_index],
    ripqp_mono1.dual_feas[pb_index],
  ]
  data2[4 * (pb_index-1) + 4, :] .= [
    "multik nops",
    ripqp_multi_nops3.elapsed_time[pb_index],
    ripqp_multi_nops3.iters_sp[pb_index],
    ripqp_multi_nops3.iters_sp2[pb_index],
    ripqp_multi_nops3.objective[pb_index],
    ripqp_multi_nops3.pdd[pb_index],
    ripqp_multi_nops3.primal_feas[pb_index],
    ripqp_multi_nops3.dual_feas[pb_index],
  ]
end

table = pretty_table(data2; 
    header = header2,
    row_names= row_names2,
    title = "bm ripqp quad prec",
    body_hlines = [4, 8],
    backend = Val(:latex),
    formatters = (ft_printf(["%7.1e", "%d", "%d", "%7.1e","%7.1e","%7.1e","%7.1e"], 2:8),
      (v, i, j) -> (SolverBenchmark.safe_latex_AbstractFloat(v)),
      )
  )