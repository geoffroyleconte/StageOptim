using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, BenchmarkProfiles, Plots
using DelimitedFiles, CSV

# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using RipQP

# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
path_pb = "/home/gelecd/.julia/artifacts/545f8c5577a056981a21caf3f53bd7b59cf67410/optrove-netlib-lp-f83996fca937"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
# save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks"
# save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/perf_profiles/lp2"
# save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/amdahl_benchmarks/perf_profiles/test_lp2"
save_path = "/home/gelecd/code/docGL/benchmarks/frontal22_results/prof2_lp"
# save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\docGL\amdahl_benchmarks\perf_profiles\test3"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
# qm = QuadraticModel(readqps(string(path_pb, "\\irish-electricity.mps")))

function createQuadraticModel(qpdata; name="qp_pb")
  return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
          Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
          lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
          c0=qpdata.c0, name=name)
end

function zeros_logscale!(v, min_val)
  for i=1:length(v)
    if v[i] == 0
      v[i] += min_val
    end
  end
end

function get_sp(formulation::Symbol, kmethod::Symbol, preconditioner)
  if kmethod ∈ [:tricg, :trimr, :gpmr, :lslq, :lsqr, :lsmr, :lnlq, :craig, :craigmr]
    sp = eval(formulation)(kmethod = kmethod, atol0 = 0.1, rtol0 = 0.1,
                          atol_min=1.0e-3, rtol_min=1.0e-1,
                          ρ_min = 1e1 * sqrt(eps()), δ_min = 1e1 * sqrt(eps()),
                          mem = 100, k3_resid = true, cb_only = true,
                          )
  else
    sp = eval(formulation)(kmethod = kmethod, preconditioner = preconditioner,
                          atol0 = 0.1, rtol0 = 0.1,
                          atol_min=1.0e-3, rtol_min=1.0e-1,
                          ρ_min = 1e1 * sqrt(eps()), δ_min = 1e1 * sqrt(eps()),
                          mem = 100, k3_resid = true, cb_only = true,
                          )
  end
  return sp
end

function ripqp_generic_solver(qm, sp)
  return RipQP.ripqp(qm, display = false, sp = sp, solve_method = IPF(), 
                     #, stepsize = stepsize,
                     # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                     itol = RipQP.InputTol(max_iter=200, max_time=30.0,
                                           ϵ_pdd = 1.0e-4, ϵ_rb = 1.0e-4, ϵ_rc = 1.0e-4,
                                           ))
end

# init solvers list
sps = []
push!(sps, K2LDLParams())
solvers = [:K2_LDL]
global c_solv = 1
# K1
for kmethod ∈ [:minres, :cg, :cg_lanczos, :cr, :bicgstab, :qmr, :diom, :dqgmres, :gmres] #:fom,
  push!(sps, get_sp(:K1KrylovParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K1_, kmethod))
end

# # K1.1 Structured
# for kmethod ∈ [:lslq, :lsqr, :lsmr]
#   push!(sps, get_sp(:K1_1StructuredParams, kmethod, :Identity))
#   global c_solv += 1
#   push!(solvers, Symbol(:K1_1_, kmethod))
# end

# # K1.2 Structured
# for kmethod ∈ [:lnlq, :craig, :craigmr]
#   push!(sps, get_sp(:K1_2StructuredParams, kmethod, :Identity))
#   global c_solv += 1
#   push!(solvers, Symbol(:K1_2_, kmethod))
# end

# K2
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres] #:fom,
  push!(sps, get_sp(:K2KrylovParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K2_, kmethod))
end

# K2 Jacobi
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres]
  push!(sps, get_sp(:K2KrylovParams, kmethod, Jacobi()))
  global c_solv += 1
  push!(solvers, Symbol(:K2Jacobi_, kmethod))
end

# K2 Equilibration
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres]
  push!(sps, get_sp(:K2KrylovParams, kmethod, Equilibration()))
  global c_solv += 1
  push!(solvers, Symbol(:K2Equilibration_, kmethod))
end

# # K2 structured
for kmethod ∈ [:tricg, :trimr] #, :gpmr
  push!(sps, get_sp(:K2StructuredParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K2_, kmethod))
end

# # K2.5
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres] #:fom, 
  push!(sps, get_sp(:K2_5KrylovParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K2_5_, kmethod))
end

# # K2.5 Jacobi
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres] #:fom, 
  push!(sps, get_sp(:K2_5KrylovParams, kmethod, Jacobi()))
  global c_solv += 1
  push!(solvers, Symbol(:K2_5Jacobi_, kmethod))
end

# # K2.5 Structured
for kmethod ∈ [:tricg, :trimr] #, :gpmr
  push!(sps, get_sp(:K2_5StructuredParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K2_5_, kmethod))
end

# # K3
for kmethod ∈ [:bicgstab, :qmr, :diom, :dqgmres, :gmres] # :fom, 
  push!(sps, get_sp(:K3KrylovParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K3_, kmethod))
end

# # K3S
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres] #:fom, 
  push!(sps, get_sp(:K3SKrylovParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K3S_, kmethod))
end

# # K3S Structured
# for kmethod ∈ [:tricg, :trimr, :gpmr]
#   push!(sps, get_sp(:K3SStructuredParams, kmethod, :Identity))
#   global c_solv += 1
#   push!(solvers, Symbol(:K3S_, kmethod))
# end

# # K3.5
for kmethod ∈ [:minres, :minres_qlp, :symmlq, :diom, :dqgmres, :gmres]# :fom,
  push!(sps, get_sp(:K3_5KrylovParams, kmethod, Identity()))
  global c_solv += 1
  push!(solvers, Symbol(:K3_5_, kmethod))
end

# # K3.5 Structured
# for kmethod ∈ [:tricg, :trimr, :gpmr]
#   push!(sps, get_sp(:K3_5StructuredParams, kmethod, :Identity))
#   global c_solv += 1
#   push!(solvers, Symbol(:K3_5_, kmethod))
# end

function get_QuadraticModels(path_pb :: String, n_pb :: Int)
  qms = []
  i_max = n_pb
  i = 1
  for file_name in readdir(path_pb)
    if file_name[end-3:end] == ".SIF"
      println(file_name)
      pb_i = string(path_pb, "/", file_name)
      if file_name in ["BLEND.SIF"; "DFL001.SIF"; "FORPLAN.SIF"; "GFRD-PNC.SIF"; "SIERRA.SIF";
                  "EXDATA.SIF"; "QFORPLAN.SIF"; "QGFRDXPN.SIF"; "VALUES.SIF"]
          qpdata_i = readqps(pb_i, mpsformat=:fixed)
      else
          qpdata_i = readqps(pb_i)
      end
      qm = createQuadraticModel(qpdata_i, name=file_name[1:end-4])
      push!(qms, qm)
      if i==i_max
        break
      end
      i += 1
    end
  end
  return qms
end

function optimize_ripqp!(qms, ripqp_func :: Function, data, solver)
  n_pb = length(qms)
  n_iter = size(data)[1]
  for i=1:n_pb
    qm = qms[i]
    stats = ripqp_func(qm)
    ldata = length(stats.solver_specific[:pddH])
    data[1:ldata, i, solver] .= stats.solver_specific[:pddH]
    if ldata < n_iter
      data[ldata+1:end, i, solver] .= stats.solver_specific[:pddH][end]
    end
  end
end

n_pb = 1000

pb_i = string(path_pb, "/", "AFIRO.SIF")
qpdata_i = readqps(pb_i)
qm = createQuadraticModel(qpdata_i, name="AFIRO")

fsolvers = Vector{Function}(undef, c_solv)
for i=1:length(sps)
  spi = sps[i]
  solver_i(qm) = ripqp_generic_solver(qm, spi)
  fsolvers[i] = solver_i
end

for i=1:length(fsolvers)
  solver_i = fsolvers[i]
  stats_compile = solver_i(qm)
  println(string(solvers[i]), " compiled")
end

n_solvers = length(solvers)

function save_problems(file_path :: String, ripqp_func :: Function, solv_str::String, qms)
  stats = solve_problems(ripqp_func, qms)
  # file = jldopen(string(file_path, "/", solv_str, ".jld2"), "w")
  # file["stats"] = stats
  # close(file)
  CSV.write(string(file_path, "/", solv_str, ".csv"), stats)
  return Nothing
end

qms = get_QuadraticModels(path_pb, n_pb)
solvers_list = [string(solver) for solver in solvers]
for is in 1: length(fsolvers)
  save_problems(save_path, fsolvers[is], solvers_list[is], qms)
  println(string(solvers[is]), " done")
end

open(string(save_path, "/test2_solvs.txt"), "w") do io
  writedlm(io, solvers_list)
end;