using RegularizedOptimization, RegularizedProblems, ProximalOperators, NLPModels, NLPModelsModifiers, LinearAlgebra

function unbounded_model(β, α, ϵ, γ3, p)

  @assert β ≥ 2 / α + 1
  kϵ = Int(floor(ϵ^(-2 / (1 - p))))
  K = 0:kϵ # K[k] = k - 1
  w(k) = (kϵ - k) / kϵ 
  Δ(k) = γ3^k
  G = [- ϵ * (1 + w(k)) for k in K] # G[k] = g_(k-1)
  g0 = -2 * ϵ
  G[1] = g0
  S = [-1 / k^p * G[k+1] for k in K] # S[k] = s_(k-1)
  S[1] = -g0
  X = similar(S) # X[k] = x_(k-1)
  X[1] = 0.0
  for k=1:kϵ
    X[k+1] = X[k] + S[k]
  end
  f0 = 8 * ϵ^2 + 4.0 / (1.0 - p)
  F = zeros(kϵ + 1) # F[k] = f_(k-1)
  F[1] = f0
  for k=1:kϵ
    F[k+1] = F[k] + G[k] * S[k] 
  end
  s_m1 = 100 * abs(g0)

  A = [1.0 1.0;
      2.0 3.0]
  As = copy(A)
  rhs = zeros(2)
  # π (τ) = ∑²ᵖ⁺¹ cᵢ τⁱ
  # cᵢ = f₀⁽ⁱ⁾ / i!  pour i ∈ [ 0 : p]
  C = [zeros(4) for k in K] # C[k] = coeffs π_(k-1) sur [x_(k-1); x_k]
  for k=1:kϵ
    C[k][1] = F[k]
    C[k][2] = G[k]
    sk = S[k]
    As[1, 1] = A[1, 1] * sk^2
    As[1, 2] = A[1, 2] * sk^3
    As[2, 1] = A[2, 1] * sk
    As[2, 2] = A[2, 2] * sk^2
    rhs[2] = G[k+1] - G[k]
    c23 = As \ rhs
    C[k][3] = c23[1]
    C[k][4] = c23[2]
  end
  # coef [-1; 0]
  c_m1 = zeros(4)
  c_m1[1] = f0
  As[1, 1] = A[1, 1] * s_m1^2
  As[1, 2] = A[1, 2] * s_m1^3
  As[2, 1] = A[2, 1] * s_m1
  As[2, 2] = A[2, 2] * s_m1^2
  rhs[2] = G[1]
  c23 = As \ rhs
  c_m1[3] = c23[1]
  c_m1[4] = c23[2]
  # coef [kϵ; kϵ + 1]
  C[kϵ+1][1] = F[kϵ+1]
  C[kϵ+1][2] = G[kϵ+1]
  sk = S[kϵ+1]
  As[1, 1] = A[1, 1] * sk^2
  As[1, 2] = A[1, 2] * sk^3
  As[2, 1] = A[2, 1] * sk
  As[2, 2] = A[2, 2] * sk^2
  rhs[2] = - G[kϵ+1]
  c23 = As \ rhs
  C[kϵ+1][3] = c23[1]
  C[kϵ+1][4] = c23[2]

  function obj(x)
    k = findfirst(xk -> xk ≥ x[1], X) # k tel que x ∈ [X[k-1]; X[k]] soit x ∈ [x_(k-2); x_(k-1)]
    if k === nothing
      Ci = C[kϵ+1]
      τ = x[1] - X[kϵ+1]
    elseif k == 1
      if x[1] < -s_m1
        return F[1]
      else
        Ci = c_m1
        τ = x[1] + s_m1
      end
    else
      Ci = C[k-1]
      τ = x[1] - X[k-1]
    end
    return Ci[1] + Ci[2] * τ + Ci[3] * τ^2 + Ci[4] * τ^3
  end

  function grad!(g, x)
    k = findfirst(xk -> xk ≥ x[1], X)
    if k === nothing
      Ci = C[kϵ+1]
      τ = x[1] - X[kϵ+1]
    elseif k == 1
      if x[1] < -s_m1
        return 0.0
      else
        Ci = c_m1
        τ = x[1] + s_m1
      end
    else
      Ci = C[k-1]
      τ = x[1] - X[k-1]
    end
    g[1] = Ci[2] + 2 * Ci[3] * τ + 3 * Ci[4] * τ^2
    return g
  end

  return FirstOrderModel(obj, grad!, [0.0]), X, kϵ
end

β = 1e16
ϵ = 1.0e-1
p = 1 / 10
α = 1 / eps()#eps()
γ = 1.0001
ν = 1 / sqrt(1 + 2 / α)
options = ROSolverOptions(ν = ν, β = β, ϵa = ϵ, ϵr = 0.0, α = α, γ = γ, verbose = 10)
options_inner = ROSolverOptions(maxIter = 10)
model, X, kϵ = unbounded_model(β, α, ϵ, γ, p)
h = NormL1(0.0)
χ = NormLinf(1.0)
lsr1model = LSR1Model(model)
TR_out = TR(lsr1model, h, χ, options, x0 = lsr1model.meta.x0)
# TR_out = TRDH(model, h, χ, options, x0 = lsr1model.meta.x0)



γ3 = γ
@assert β ≥ 2 / α + 1
kϵ = Int(floor(ϵ^(-2 / (1 - p))))
K = 0:kϵ # K[k] = k - 1
w(k) = (kϵ - k) / kϵ 
Δ(k) = γ3^k
G = [- ϵ * (1 + w(k)) for k in K] # G[k] = g_(k-1)
g0 = -2 * ϵ
G[1] = g0
S = [-1 / k^p * G[k+1] for k in K] # S[k] = s_(k-1)
S[1] = -g0
X = similar(S) # X[k] = x_(k-1)
X[1] = 0.0
for k=1:kϵ
  X[k+1] = X[k] + S[k]
end
f0 = 8 * ϵ^2 + 4.0 / (1.0 - p)
F = zeros(kϵ + 1) # F[k] = f_(k-1)
F[1] = f0
for k=1:kϵ
  F[k+1] = F[k] + G[k] * S[k] 
end
s_m1 = 100 * abs(g0)

A = [1.0 1.0;
    2.0 3.0]
As = copy(A)
rhs = zeros(2)
# π (τ) = ∑²ᵖ⁺¹ cᵢ τⁱ
# cᵢ = f₀⁽ⁱ⁾ / i!  pour i ∈ [ 0 : p]
C = [zeros(4) for k in K] # C[k] = coeffs π_(k-1) sur [x_(k-1); x_k]
for k=1:kϵ
  C[k][1] = F[k]
  C[k][2] = G[k]
  sk = S[k]
  As[1, 1] = A[1, 1] * sk^2
  As[1, 2] = A[1, 2] * sk^3
  As[2, 1] = A[2, 1] * sk
  As[2, 2] = A[2, 2] * sk^2
  rhs[2] = G[k+1] - G[k]
  c23 = As \ rhs
  C[k][3] = c23[1]
  C[k][4] = c23[2]
end
# coef [-1; 0]
c_m1 = zeros(4)
c_m1[1] = f0
As[1, 1] = A[1, 1] * s_m1^2
As[1, 2] = A[1, 2] * s_m1^3
As[2, 1] = A[2, 1] * s_m1
As[2, 2] = A[2, 2] * s_m1^2
rhs[2] = G[1]
c23 = As \ rhs
c_m1[3] = c23[1]
c_m1[4] = c23[2]
# coef [kϵ; kϵ + 1]
C[kϵ+1][1] = F[kϵ+1]
C[kϵ+1][2] = G[kϵ+1]
sk = S[kϵ+1]
As[1, 1] = A[1, 1] * sk^2
As[1, 2] = A[1, 2] * sk^3
As[2, 1] = A[2, 1] * sk
As[2, 2] = A[2, 2] * sk^2
rhs[2] = - G[kϵ+1]
c23 = As \ rhs
C[kϵ+1][3] = c23[1]
C[kϵ+1][4] = c23[2]

function f(x)
  k = findfirst(xk -> xk ≥ x, X) # k tel que x ∈ [X[k-1]; X[k]] soit x ∈ [x_(k-2); x_(k-1)]
  if k === nothing
    Ci = C[kϵ+1]
    τ = x - X[kϵ+1]
  elseif k == 1
    if x < -s_m1
      return F[1]
    else
      Ci = c_m1
      τ = x + s_m1
    end
  else
    Ci = C[k-1]
    τ = x - X[k-1]
  end
  return Ci[1] + Ci[2] * τ + Ci[3] * τ^2 + Ci[4] * τ^3
end

function g(x)
  k = findfirst(xk -> xk ≥ x, X)
  if k === nothing
    Ci = C[kϵ+1]
    τ = x - X[kϵ+1]
  elseif k == 1
    if x[1] < -s_m1
      return 0.0
    else
      Ci = c_m1
      τ = x + s_m1
    end
  else
    Ci = C[k-1]
    τ = x - X[k-1]
  end
  return Ci[2] + 2 * Ci[3] * τ + 3 * Ci[4] * τ^2
end


using PGFPlots
absc = X[1]:0.01:X[end]
b = Axis(
  PGFPlots.Plots.Linear(absc, f.(absc), mark = "none"),
  xlabel = "x",
  ylabel = "f(x)",
  # ymode = "log",
)
# save("f-plot.tikz", b)

c = Axis(
  PGFPlots.Plots.Linear(absc, g.(absc), mark = "none"),
  xlabel = "x",
  ylabel = "f'(x)",
  # ymode = "log",
)
# save("g-plot.tikz", c)