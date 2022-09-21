res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\ripqp_paper"
using Plots
using DataFrames, SolverBenchmark, SolverTools
using JLD2
using CSV
# using FileIO
open_file(fname; res_path = res_path) = CSV.read(string(res_path, "\\", fname, ".csv"), DataFrame)

ripqp_multik1 = open_file("ripqp_multik1_quad") # 3 sp init 
ripqp_multi1 = open_file("ripqp_multi1_quad")