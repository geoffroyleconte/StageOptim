include("data_def.jl")
using LinearOperators, QuadraticModels, SparseArrays

# https://jamesmccaffrey.wordpress.com/2018/03/14/datasets-for-binary-classification/
# https://machinelearningmastery.com/standard-machine-learning-datasets/

function mulopK!(res, v, α, β::T, Kfunc::Function, X, y) where T
  n = length(res)
  res .= β == zero(T) ? zero(T) : β .* res 
  for i=1:n
    for j=i+1:n
      Kij = @views Kfunc(X[i, :], X[j, :])
      res[i] += α * Kij * y[i] * v[j]
      res[j] += α * Kij * y[j] * v[i]
    end
  end
end

K_train = zeros(n_dat_train, n_dat_train)
for i=1:n_dat_train
  for j=1:i
    K_train[i, j] = @views Krbf(X_train[i,:], X_train[j,:])
  end
end

function mulopK2!(res::AbstractVector{T}, v, K, y) where T
  mul!(res, K, v .* y)
  res .*= y
end

opK2 = LinearOperator(Float64, nvar, nvar, true, true, (res, v) -> mulopK2!(res, v, K_train, y_train))

opK = LinearOperator(Float64, nvar, nvar, true, true, (res, v, α, β) -> mulopK!(res, v, α, β, Krbf, X_train, y_train))

include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")

qm = QuadraticModel(c, Symmetric(tril!(Diagonal(y_train) * K_train * Diagonal(y_train)), :L), A=A, lcon = b, ucon = b, lvar = lvar, uvar = uvar)
# qm = QuadraticModel(c, sparse(tril!(Diagonal(y_train) * K_train * Diagonal(y_train))), A=A, lcon = b, ucon = b, lvar = lvar, uvar = uvar)
# qm = QuadraticModel(c, opK2, A=A, lcon = b, ucon = b, lvar = lvar, uvar = uvar)

stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                        # sp = RipQP.K2LDLParams(),
                        sp = RipQP.K2LDLDenseParams(fact_alg = :bunchkaufman),
                        solve_method=:IPF, scaling = false, history=false, presolve=false,
                        # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                        ),
                     itol = RipQP.InputTol(max_iter=50, max_time=20.0,
                     ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                     ))
println(stats1)

n_dat_test = length(y_test)
K_pred = zeros(n_dat_train)
function update_Kpred!(K_pred::AbstractVector{T}, Kfunc::Function, x_pred::AbstractVector{T}, X_train::AbstractMatrix{T}) where T
  for i=1:length(K_pred)
    K_pred[i] = @views Kfunc(x_pred, X_train[i, :])
  end
end

function predict_sample(α::AbstractVector{T}, K_pred::AbstractVector{T}, Kfunc::Function, 
                        x_pred::AbstractVector{T}, X_train::AbstractMatrix{T}, y_train::AbstractVector, b::T) where T
  update_Kpred!(K_pred, Kfunc, x_pred, X_train)
  K_pred .*= y_train .* α
  return sign(sum(K_pred) + b) 
end

function predict_ripqp(α::AbstractVector{T}, Kfunc, X_train, y_train, X_test, b) where {T}
  n_dat_test = size(X_test, 1)
  y_pred = zeros(n_dat_test)
  K_pred = zeros(n_dat_train)
  for i=1:n_dat_test
    y_pred[i] = @views predict_sample(α, K_pred, Kfunc, X_test[i, :], X_train, y_train, b)
  end
  return y_pred
end

y_pred1 =  predict_ripqp(stats1.solution, Krbf, X_train, y_train, X_test, -stats1.multipliers[1])
println("accuracy ripqp = " , sum(y_pred1 .== y_test) / n_dat_test)
