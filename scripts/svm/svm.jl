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
  for j=1:n_dat_train
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

qm = QuadraticModel(c, sparse(tril!(Diagonal(y_train) * K_train * Diagonal(y_train))), A=A, lcon = b, ucon = b, lvar = lvar, uvar = uvar)

# qm = QuadraticModel(c, opK2, A=A, lcon = b, ucon = b, lvar = lvar, uvar = uvar)

stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                        # sp = RipQP.K2KrylovParams(kmethod=:minres,
                        # atol0 = 0.1, rtol0 = 0.1, atol_min=1.0e-10, rtol_min=1.0e-10), 
                        solve_method=:IPF, scaling = true, history=false, presolve=false,
                        # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                        ),
                     itol = RipQP.InputTol(max_iter=50, max_time=20.0,
                     ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                     ))
println(stats1)

function predict_ripqp(α::AbstractVector{T}, x_pred, X, y, K_func) where {T}
  sumpred = b_offset
  for i=1:length(α)
    sumpred += @views α[i] * y[i] * K_func(x_pred, X[i, :])
  end
  return sign(sumpred)
end
