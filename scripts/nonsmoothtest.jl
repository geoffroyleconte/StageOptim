using Random
using LinearAlgebra
using ProximalOperators
using NLPModels, NLPModelsModifiers, RegularizedProblems, RegularizedOptimization, ShiftedProximalOperators
using Printf

# compound = 1
# model, nls_model, sol = bpdn_model(compound)
# f = LSR1Model(model)
# λ = norm(grad(model, zeros(model.meta.nvar)), Inf) / 10
# h = NormL1(λ)
# options = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-6, ϵr = 1e-6, verbose = 10, maxIter = 100)
# # χ = NormL2(1.0)
# χ = NormLinf(1.0)
# x = rand(10)
# # shifted(h, x, 0.5, χ)
# reset!(f)
# TR_out = TR(f, h, χ, options, x0 = f.meta.x0)
# # R2_out = R2(f, h, options, x0 = f.meta.x0)
# # TRDH_out = TRDH(f, h, χ, options, x0 = f.meta.x0, spectral = false, psb = true)

# # nls
# model, nls_model, sol = bpdn_model(compound, bounds = true)
# λ = norm(grad(model, zeros(model.meta.nvar)), Inf) / 10
# h = NormL1(λ)
# options = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-6, ϵr = 1e-6, verbose = 10, maxIter = 100)
# sub_options = ROSolverOptions(spectral = true)
# # χ = NormL2(1.0)
# χ = NormLinf(1.0)
# x = rand(10)
# reset!(nls_model)
# LMTR_out = LMTR(nls_model, h, χ, options, x0 = nls_model.meta.x0, subsolver = R2, subsolver_options = sub_options)

# cstr
# Random.seed!(1010)
# compound = 1
# model, nls_model, sol = bpdn_model(compound, bounds = true)
# f = LSR1Model(model)
# λ = norm(grad(model, zeros(model.meta.nvar)), Inf) / 2
# h = NormL1(λ)
# χ = NormLinf(1.0)
# opt_RIPM = RIPMOptions(
#   maxIter_inner = 40,
#   μ0 = 1.0e-3,
#   μmin = sqrt(eps()) * 10,
#   ϵ0 = 1.0,
#   δ0 = 1.0,
#   δmin = sqrt(eps()),
#   resetQN = true,
# )
# options = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-4, ϵr = 1e-4, verbose = 10, maxIter = 40, opt_RIPM = opt_RIPM)
# subsolver_options = ROSolverOptions(maxIter = 100)
# @info " using TR to solve with" h χ
# reset!(f)
# TR_out = TR(f, h, χ, options, x0 = f.meta.x0)
# RIPM_out = RIPM(f, h, χ, options, x0 = f.meta.x0, subsolver_options = subsolver_options)
# R2_out = R2(f, h, options, x0 = f.meta.x0)
# TRDH_out = TRDH(f, h, χ, options, x0 = f.meta.x0, spectral = true, psb = false)

## FH cstr
# using ADNLPModels, DifferentialEquations
# Random.seed!(1234)
# cstr = true
# lvar = [-Inf, 0.5, -Inf, -Inf, -Inf]
# uvar = fill(Inf, 5)
# data, simulate, resid, misfit, x0 = RegularizedProblems.FH_smooth_term()
# model = ADNLPModel(misfit, ones(5), lvar, uvar)
# f = LBFGSModel(model)
# λ = 2.0e1
# h = NormL1(λ)
# χ = NormLinf(1.0)
# opt_RIPM = RIPMOptions(
#   maxIter_inner = 40,
#   μ0 = 1.0e-3,
#   μmin = sqrt(eps()) * 10,
#   ϵ0 = 1.0,
#   δ0 = 1.0,
#   δmin = sqrt(eps()),
#   resetQN = true,
# )
# options = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-4, ϵr = 1e-4, verbose = 10, maxIter = 100, opt_RIPM = opt_RIPM)
# subsolver_options = ROSolverOptions(maxIter = 500)
# @info " using RIPM to solve with" h χ
# reset!(f)
# RIPM_out = RIPM(f, h, χ, options, x0 = f.meta.x0, subsolver_options = subsolver_options)
# reset!(f)
# TR_out = TR(f, h, χ, options, x0 = f.meta.x0, subsolver_options = subsolver_options)

## NNMF
m, n, k = 100, 50, 5
Random.seed!(1234)
model, A, selected = nnmf_model(m, n, k)
λ = 1.0e-1
h = NormL1(λ)
f = LSR1Model(model)
χ = NormLinf(1.0)
opt_RIPM = RIPMOptions(
  maxIter_inner = 40,
  μ0 = 1.0e-3,
  μmin = sqrt(eps()) * 10,
  ϵ0 = 1.0,
  δ0 = 1.0,
  δmin = sqrt(eps()),
  resetQN = true,
)
options = ROSolverOptions(ν = 1.0, β = 1e16, ϵa = 1e-4, ϵr = 1e-4, verbose = 10, maxIter = 100, opt_RIPM = opt_RIPM)
subsolver_options = ROSolverOptions(maxIter = 500)
@info " using RIPM to solve with" h χ
reset!(f)
RIPM_out = RIPM(f, h, χ, options, x0 = f.meta.x0, subsolver_options = subsolver_options, selected = selected)
reset!(f)
TR_out = TR(f, h, χ, options, x0 = f.meta.x0, subsolver_options = subsolver_options, selected = selected)