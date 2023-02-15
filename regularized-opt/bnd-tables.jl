using PrettyTables
using Random
using LinearAlgebra
using ProximalOperators
using NLPModels, NLPModelsModifiers, RegularizedProblems, RegularizedOptimization, ShiftedProximalOperators, SolverBenchmark
using Printf

Random.seed!(1234)
compound = 1

# utils for extracting stats / display table
modelname(nlp::LSR1Model) = "LSR1"
subsolvername(subsolver::Symbol) = subsolver == :None ? "" : string("-", subsolver)
function options_str(options::ROSolverOptions, solver::Symbol, subsolver_options::ROSolverOptions, subsolver::Symbol)
  if solver == :TRDH
    return !options.spectral ? (options.psb ? "-DiagQN-PSB" : "-DiagQN-Andrei") : "-Spectral"
  elseif solver == :TR && subsolver == :TRDH
    return !subsolver_options.spectral ? (subsolver_options.psb ? "-DiagQN-PSB" : "-DiagQN-Andrei") : "-Spectral"
  else
    return ""
  end
end
grad_evals(nlp::AbstractNLPModel) = neval_grad(nlp)
grad_evals(nls::AbstractNLSModel) = neval_jtprod_residual(nls) + neval_jprod_residual(nls)
function nb_prox_evals(stats, solver::Symbol)
  if solver ∈ [:TR, :R2, :TRDH]
    prox_evals = sum(stats.solver_specific[:SubsolverCounter])
  else
    error("not implemented")
  end
  return prox_evals
end

function benchmark_table(
  f::AbstractNLPModel,
  selected,
  sol,
  h,
  λ,
  solvers,
  subsolvers,
  solver_options,
  subsolver_options,
  pb_name::String,
)

  row_names = ["$(solver)$(subsolvername(subsolver))$(options_str(opt, solver, subsolver_opt, subsolver))" 
                for (solver, opt, subsolver, subsolver_opt) in zip(solvers, solver_options, subsolvers, subsolver_options)]

  n∇f_evals = []
  nprox_evals = []
  solver_stats = []

  for (solver, subsolver, opt, sub_opt) in zip(solvers, subsolvers, solver_options, subsolver_options)
    @info " using $solver with subsolver = $subsolver"
    args = solver == :R2 ? () : (NormLinf(1.0),)
    if subsolver == :None
      solver_out = eval(solver)(f, h, args..., opt, x0 = f.meta.x0, selected = selected)
    else
      solver_out = eval(solver)(f, h, args..., opt, x0 = f.meta.x0, subsolver = eval(subsolver),
                                subsolver_options = sub_opt, selected = selected)
    end
    push!(n∇f_evals, grad_evals(f))
    push!(nprox_evals, nb_prox_evals(solver_out, solver))
    push!(solver_stats, solver_out)
    reset!(f)

  end

  if length(sol) == 0
    header = [
      "f(x)",
      "h(x)/λ",
      "∇f evals",
      "prox calls",
    ]
  else
    header = [
      "f(x) (true = $(round(obj(model, sol); sigdigits = 4)))",
      "h(x)/λ",
      "||x-x*||/||x*||",
      "∇f evals",
      "prox calls",
    ]
  end

  n_solvers = length(row_names)
  data = Matrix{Any}(undef, n_solvers, length(header))
  for i=1:n_solvers
    solver_out = solver_stats[i]
    x = solver_out.solution
    fx = solver_out.solver_specific[:Fhist][end]
    hx = solver_out.solver_specific[:Hhist][end]
    n∇f = n∇f_evals[i]
    nprox = nprox_evals[i]
    if length(sol) == 0
      data[i, :] .= [fx, hx / λ, n∇f, nprox]
    else
      err = norm(x - sol) / norm(sol)
      data[i, :] .= [fx, hx / λ, err, n∇f, nprox]
    end
  end

  if length(sol) == 0
    print_formats = ft_printf(["%7.3e", "%7.1e", "%i","%i"], 1:length(header))
  else
    print_formats = ft_printf(["%7.3e", "%7.1e", "%7.3e", "%i","%i"], 1:length(header))
  end

  return pretty_table(data; 
      header = header,
      row_names= row_names,
      title = "$pb_name $(modelname(f)) $(typeof(h).name.name)",
      # backend = Val(:latex),
      formatters = (print_formats,
        # (v, i, j) -> (SolverBenchmark.safe_latex_AbstractFloat(v)),
        ),
    )
end


# model
model, nls_model, sol = bpdn_model(compound, bounds = true)

# parameters
solvers = [:R2, :TRDH, :TRDH, :TR, :TR, :TR, :TR]
subsolvers = [:None, :None, :None, :R2, :TRDH, :TRDH, :TRDH]
f = LSR1Model(model)
λ = norm(grad(model, zeros(model.meta.nvar)), Inf) / 10
h = NormL0(λ)

options = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-6, ϵr = 1e-6, verbose = 10, maxIter = 500, spectral = true)
options2 = ROSolverOptions(spectral = false, psb = true)
options3 = ROSolverOptions(spectral = false, psb = false)
options4 = ROSolverOptions(spectral = true)
options5 = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-6, ϵr = 1e-6, verbose = 10, maxIter = 500, spectral = false, psb = true)
solver_options = [options, options, options5, options, options, options, options]
subsolver_options = [options2, options2, options2, options2, options2, options3, options4] # n'importe lequel si subsolver = :None

benchmark_table(f, 1:f.meta.nvar, sol, h, λ, solvers, subsolvers, solver_options, subsolver_options, "BPDN-cstr")


# no cstr model

# model
model, nls_model, sol = bpdn_model(compound, bounds = false)
f = LSR1Model(model)
benchmark_table(f, 1:f.meta.nvar, sol, h, λ, solvers, subsolvers, solver_options, subsolver_options, "BPDN")

# NNMF
m, n, k = 100, 50, 5
model, A, selected = nnmf_model(m, n, k)
f = LSR1Model(model)
λ = norm(grad(model, rand(model.meta.nvar)), Inf) / 200
h = NormL0(λ)
options = ROSolverOptions(ν = 1.0e-3, β = 1e16, ϵa = 1e-6, ϵr = 1e-6, verbose = 10, maxIter = 500, spectral = true)
options5 = ROSolverOptions(ν = 1.0e-3, β = 1e16, ϵa = 1e-6, ϵr = 1e-6, verbose = 10, maxIter = 500, spectral = false)
solver_options = [options, options, options5, options, options, options, options]
subsolver_options = [options2, options2, options2, options2, options2, options3, options4]
benchmark_table(f, selected, [], h, λ, solvers, subsolvers, solver_options, subsolver_options,
                "NNMF with m = $m, n = $n, k = $k, ν = 1.0e-3,")







# header = ["TR LSR1 L0Box", "R2 LSR1 L0Box", "LM L0Box", "LMTR L0Box"]
# TR_out = TR(f, h, χ, options, x0 = f.meta.x0)
# n∇f_TR = neval_grad(f)
# prox_evals_TR = sum(TR_out.solver_specific[:SubsolverCounter])
# reset!(f)
# R2_out = R2(f, h, options, x0 = f.meta.x0)
# n∇f_R2 = neval_grad(f)
# prox_evals_R2 = R2_out.iter
# reset!(f)
# LM_out = LM(nls_model, h, options, x0 = nls_model.meta.x0)
# n∇f_LM = neval_jtprod_residual(nls_model) + neval_jprod_residual(nls_model)
# prox_evals_LM = sum(LM_out.solver_specific[:SubsolverCounter])
# reset!(nls_model)
# LMTR_out = LMTR(nls_model, h, χ, options, x0 = nls_model.meta.x0)
# n∇f_LMTR = neval_jtprod_residual(nls_model) + neval_jprod_residual(nls_model)
# prox_evals_LMTR = sum(LMTR_out.solver_specific[:SubsolverCounter])
# reset!(nls_model)
# n∇f_evals = [n∇f_TR, n∇f_R2, n∇f_LM, n∇f_LMTR]
# nprox_evals = [prox_evals_TR, prox_evals_R2, prox_evals_LM, prox_evals_LMTR]

# solver_stats = [TR_out, R2_out, LM_out, LMTR_out]