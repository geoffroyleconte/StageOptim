using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, Plots, Printf
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
# using RipQP
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
# qm = QuadraticModel(readqps(string(path_pb, "\\irish-electricity.mps")))
problem = "\\KB2.SIF"

function zeros_logscale!(v, min_val)
  for i=1:length(v)
    if v[i] == 0
      v[i] += min_val
    end
  end
end

max_iter, min_val, formul = 50, 1.0e-15, :IPF
kmethod, precond = :minres, :Identity
prc = plot(yaxis=:log10)
title!(prc, string("rc residuals ", problem[2:end-4]))
prb = plot(yaxis=:log10)
title!(prb, string("rb residuals ", problem[2:end-4]))
ppdd = plot(yaxis=:log10)
title!(ppdd, string("pdd residuals ", problem[2:end-4]))
pnp = plot()
title!(pnp, string("number of products ", problem[2:end-4]))
pμ = plot(yaxis=:log10)
title!(pμ, string("μ values ", problem[2:end-4]))
pmbd = plot(yaxis=:log10)
title!(pmbd, string("min distance to bounds ", problem[2:end-4]))
pkrN = plot(yaxis=(:log10))
title!(pkrN, string("norm of residuals Krylov method ", problem[2:end-4]))
pkrPN = plot(yaxis=(:log10))
title!(pkrPN, string("norm of primal residuals Krylov method ", problem[2:end-4]))
pkrDN = plot(yaxis=(:log10))
title!(pkrDN, string("norm of dual residuals Krylov method ", problem[2:end-4]))

for (ρ_min, δ_min) in [(1e2 * sqrt(eps()), 1e3 * sqrt(eps())) , (1e-5 * sqrt(eps()), 1e-5 * sqrt(eps()))]
  sρ = @sprintf("%.0e", ρ_min)
  sδ = @sprintf("%.0e", δ_min)
  qm = QuadraticModel(readqps(string(path_pb, problem), mpsformat=:fixed));
  stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                          sp = RipQP.K2KrylovParams(kmethod=kmethod, preconditioner = precond, 
                          ρ_min = ρ_min, δ_min = δ_min,
                          atol_min=1.0e-10, rtol_min=1.0e-10), 
                          solve_method=formul, history=true
                          # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                          ),
                      itol = RipQP.InputTol(max_iter=max_iter, max_time=20.0,
                      ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                      ))

  plot!(prc, stats1.solver_specific[:rcNormH], label = string("rcNorm K2 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(prb, stats1.solver_specific[:rbNormH], label = string("rbNorm K2 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(ppdd, stats1.solver_specific[:pddH], label = string("pdd K2 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(pnp, stats1.solver_specific[:nprodH], label = string("nprod K2 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(pμ, stats1.solver_specific[:μH], label = string("μ K2 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(pmbd, stats1.solver_specific[:min_bound_distH], label = string("mbd K2 ρ_min = ", sρ, " δ_min = ", sδ))
  zeros_logscale!(stats1.solver_specific[:KresNormH], min_val)
  plot!(pkrN, stats1.solver_specific[:KresNormH], label = string("krN K2 ρ_min = ", sρ, " δ_min = ", sδ))

  qm = QuadraticModel(readqps(string(path_pb, problem), mpsformat=:fixed));
  stats2 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                          sp = RipQP.K2_5KrylovParams(kmethod=kmethod, preconditioner = precond,
                          ρ_min = ρ_min, δ_min = δ_min,
                          atol_min=1.0e-10, rtol_min=1.0e-10), 
                          solve_method=formul, history=true
                          # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                          ),
                      itol = RipQP.InputTol(max_iter=max_iter, max_time=20.0,
                      ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                      ))

  plot!(prc, stats2.solver_specific[:rcNormH], label = string("rcNorm K2.5 ρ_min = ", sρ , " δ_min = ", sδ))
  plot!(prb, stats2.solver_specific[:rbNormH], label = string("rbNorm K2.5 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(ppdd, stats2.solver_specific[:pddH], label = string("pdd K2.5 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(pnp, stats2.solver_specific[:nprodH], label = string("nprod K2.5 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(pμ, stats2.solver_specific[:μH], label = string("μ K2.5 ρ_min = ", sρ, " δ_min = ", sδ))
  plot!(pmbd, stats2.solver_specific[:min_bound_distH], label = string("mbd K2.5 ρ_min = ", sρ, " δ_min = ", sδ))
  zeros_logscale!(stats2.solver_specific[:KresNormH], min_val)
  plot!(pkrN, stats2.solver_specific[:KresNormH], label = string("krN K2.5 ρ_min = ", sρ, " δ_min = ", sδ))
end
plot!(ppdd)
