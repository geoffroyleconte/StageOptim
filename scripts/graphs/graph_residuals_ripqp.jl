using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, Plots
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
# using RipQP
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
# qm = QuadraticModel(readqps(string(path_pb, "\\irish-electricity.mps")))
problem = "\\AGG.SIF"

function zeros_logscale!(v, min_val)
  for i=1:length(v)
    if v[i] == 0
      v[i] += min_val
    end
  end
end

max_iter, min_val, formul = 10, 1.0e-15, :IPF
condshift = 1
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
pAcond = plot(yaxis=(:log10))
title!(pAcond, string("Acond Krylov method ", problem[2:end-4]))

for kmethod in [:minres]#, :minres_qlp]#, :minares]
  for precond in [:Identity, :Jacobi]
    qm = QuadraticModel(readqps(string(path_pb, problem), mpsformat=:fixed));
    stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                            sp = RipQP.K2KrylovParams(kmethod=kmethod, preconditioner = precond, atol_min=1.0e-10, rtol_min=1.0e-10), 
                            solve_method=formul, history=true
                            # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                            ),
                        itol = RipQP.InputTol(max_iter=max_iter, max_time=20.0,
                        ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                        ))

    plot!(prc, stats1.solver_specific[:rcNormH], label = string("rcNorm K2 ", kmethod, " ", precond))
    plot!(prb, stats1.solver_specific[:rbNormH], label = string("rbNorm K2 ", kmethod, " ", precond))
    plot!(ppdd, stats1.solver_specific[:pddH], label = string("pdd K2 ", kmethod, " ", precond))
    plot!(pnp, stats1.solver_specific[:nprodH], label = string("nprod K2 ", kmethod, " ", precond))
    plot!(pμ, stats1.solver_specific[:μH], label = string("μ K2 ", kmethod, " ", precond))
    plot!(pmbd, stats1.solver_specific[:min_bound_distH], label = string("mbd K2 ", kmethod, " ", precond))
    zeros_logscale!(stats1.solver_specific[:KresNormH], min_val)
    plot!(pkrN, stats1.solver_specific[:KresNormH], label = string("krN K2 ", kmethod, " ", precond))
    plot!(pAcond, stats1.solver_specific[:AcondH] .+condshift, label = string("pAcond K2 ", kmethod, " ", precond))

    qm = QuadraticModel(readqps(string(path_pb, problem), mpsformat=:fixed));
    str2 = "K2.5 "
    stats2 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                            sp = RipQP.K2_5KrylovParams(kmethod=kmethod, preconditioner = precond, 
                            # atol0 = 0., rtol0 = 0.1,
                            atol_min=1.0e-10, rtol_min=1.0e-10), 
                            solve_method=formul, history=true,
                            # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                            ),
                        itol = RipQP.InputTol(max_iter=max_iter, max_time=20.0,
                        ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                        ))

    plot!(prc, stats2.solver_specific[:rcNormH], label = string("rcNorm ", str2, kmethod, " ", precond))
    plot!(prb, stats2.solver_specific[:rbNormH], label = string("rbNorm ", str2, kmethod, " ", precond))
    plot!(ppdd, stats2.solver_specific[:pddH], label = string("pdd ", str2, kmethod, " ", precond))
    plot!(pnp, stats2.solver_specific[:nprodH], label = string("nprod ", str2, kmethod, " ", precond))
    plot!(pμ, stats2.solver_specific[:μH], label = string("μ ", str2, " ", precond))
    plot!(pmbd, stats2.solver_specific[:min_bound_distH], label = string("mbd ", str2, kmethod, " ", precond))
    zeros_logscale!(stats2.solver_specific[:KresNormH], min_val)
    plot!(pkrN, stats2.solver_specific[:KresNormH], label = string("krN ", str2, kmethod, " ", precond))
    plot!(pAcond, stats2.solver_specific[:AcondH] .+ condshift, label = string("pAcond K2.5 ", kmethod, " ", precond))
  end
end
plot!(ppdd)

# save_path = string(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\graphes\residuals","/", formul, problem[1:end-4])#, "_10-6")
# savefig(prb, string(save_path, "/prb.pdf"))
# savefig(prc, string(save_path, "/prc.pdf"))
# savefig(ppdd, string(save_path, "/ppdd.pdf"))
# savefig(pnp, string(save_path, "/pnp.pdf"))
# savefig(pμ, string(save_path, "/pmu.pdf"))
# savefig(pmbd, string(save_path, "/pmbd.pdf"))
# savefig(pkrN, string(save_path, "/pkrN.pdf"))