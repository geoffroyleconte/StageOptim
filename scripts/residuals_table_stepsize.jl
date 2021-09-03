using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, Plots, PrettyTables
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
# using RipQP

function zeros_logscale!(v, min_val)
  for i=1:length(v)
    if v[i] == 0
      v[i] += min_val
    end
  end
end

function build_table(problem)

  path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
  # path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
  save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
  # path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
  # qm = QuadraticModel(readqps(string(path_pb, "\\irish-electricity.mps")))

  stepsizes = [:classic, :Ninf]
  precond = :Identity
  max_iter = 50

  header = ["pdd", "||rb||", "||rc||", "μ", "ipm iter", "krylov iter"]
  data = zeros(length(stepsizes)*2, length(header))
  lign = 1
  row_names = []

  for kmethod in [:minres]#, :minres_qlp]#, :minares]
    for stepsize in stepsizes
      qm = QuadraticModel(readqps(string(path_pb, problem), mpsformat=:fixed));
      stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                              sp = RipQP.K2KrylovParams(kmethod=kmethod, preconditioner = precond, atol_min=1.0e-10, rtol_min=1.0e-10), 
                              solve_method=:IPF, history=true, stepsize = stepsize,
                              # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                              ),
                          itol = RipQP.InputTol(max_iter=max_iter, max_time=20.0,
                          ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                          ))

      argmin_it = argmin(stats1.solver_specific[:pddH])
      push!(row_names, "K2 stepsize $stepsize")
      push!(row_names, "K2.5 stepsize $stepsize")
      data[lign, :] = [stats1.solver_specific[:pddH][argmin_it],
                    stats1.solver_specific[:rbNormH][argmin_it],
                    stats1.solver_specific[:rcNormH][argmin_it],
                    stats1.solver_specific[:μH][argmin_it],
                    argmin_it,
                    sum(stats1.solver_specific[:nprodH][1:argmin_it])]
      lign += 1

      qm = QuadraticModel(readqps(string(path_pb, problem), mpsformat=:fixed));
      str2 = "small tol "
      stats2 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                              sp = RipQP.K2_5KrylovParams(kmethod=kmethod, preconditioner = precond, 
                                                          atol_min=1.0e-10, rtol_min=1.0e-10), 
                              solve_method=:IPF, history=true, stepsize = stepsize,
                              # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                              ),
                          itol = RipQP.InputTol(max_iter=max_iter, max_time=20.0,
                          ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                          ))

      argmin_it = argmin(stats2.solver_specific[:pddH])
          data[lign, :] = [stats2.solver_specific[:pddH][argmin_it],
                      stats2.solver_specific[:rbNormH][argmin_it],
                      stats2.solver_specific[:rcNormH][argmin_it],
                      stats2.solver_specific[:μH][argmin_it],
                      argmin_it,
                      sum(stats2.solver_specific[:nprodH][1:argmin_it])]
          lign += 1
    end
  end

  return pretty_table(data; 
                      header = header,
                      row_names= row_names,
                      title = problem[2:end-4],
                      formatters = ft_printf("%7.1e", [1, 2, 3, 4]))
end

problem = "\\AGG.SIF"
build_table(problem)

# save_path = string(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\residuals","/", formul, problem[1:end-4])#, "_10-6")