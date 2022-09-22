res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\docGL\\benchmarks\\ripqp_paper"
using Plots
using DataFrames, SolverBenchmark, SolverTools
using JLD2
using CSV
# using FileIO
open_file(fname; res_path = res_path) = CSV.read(string(res_path, "\\", fname, ".csv"), DataFrame)

ripqp_multik1 = open_file("ripqp_multik1_quad") # 3 sp init 
ripqp_multik2 = open_file("ripqp_multik2_quad") # 
ripqp_multik3 = open_file("ripqp_multik3_quad") # 1 100iter gmres ir no eq
ripqp_multik4 = open_file("ripqp_multik4_quad") # 
ripqp_multik5 = open_file("ripqp_multik5_quad") # 1 100 iter gmres no eq
ripqp_multik6 = open_file("ripqp_multik6_quad") # 
ripqp_multik7 = open_file("ripqp_multik7_quad") # 
ripqp_multi1 = open_file("ripqp_multi1_quad")