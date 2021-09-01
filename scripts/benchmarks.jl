function classic_append!(a)
    append!(a, ones(10000))
  end
  
  # BenchmarkTools.Trial: 
  # memory estimate:  1.60 MiB
  # allocs estimate:  5
  # --------------
  # minimum time:     202.900 μs (0.00% GC) 
  # median time:      305.600 μs (0.00% GC) 
  # mean time:        443.670 μs (21.98% GC)
  # maximum time:     5.279 ms (83.14% GC)  
  # --------------
  # samples:          10000
  # evals/sample:     1

function classic_append2!(a)
    T = eltype(a)
    append!(a, (-one(T) for _=1:10000))
end
 
# BenchmarkTools.Trial: 
#   memory estimate:  1.83 MiB
#   allocs estimate:  20003
#   --------------
#   minimum time:     1.390 ms (0.00% GC)
#   median time:      1.619 ms (0.00% GC)
#   mean time:        1.930 ms (5.45% GC)
#   maximum time:     33.753 ms (0.00% GC)
#   --------------
#   samples:          2590
#   evals/sample:     1
 
function new_append!(a)
    a = [a; ones(1000)]
end
 
#  BenchmarkTools.Trial: 
#    memory estimate:  1.54 MiB
#    allocs estimate:  5
#    --------------
#    minimum time:     104.500 μs (0.00% GC)
#    median time:      120.100 μs (0.00% GC)
#    mean time:        143.583 μs (9.98% GC)
#    maximum time:     3.420 ms (95.83% GC)
#    --------------
#    samples:          10000
#    evals/sample:     1
 
function new_append2!(a)
    T = eltype(a)
    n = length(a)
    a = [a; (one(T) for _=1:n)]
end
 
#   @benchmark new_append2!(ones(100000))
#  BenchmarkTools.Trial: 
#    memory estimate:  3.05 MiB
#    allocs estimate:  100025  
#    --------------
#    minimum time:     666.700 μs (0.00% GC)
#    median time:      683.799 μs (0.00% GC)
#    mean time:        873.777 μs (19.79% GC)
#    maximum time:     5.641 ms (79.81% GC)
#    --------------
#    samples:          5720
#    evals/sample:     1


using NLPModels, LinearAlgebra, Printf

function armijo(xk, dk, fk, gk, f)
    slope = dot(gk, dk) #doit être <0
    t = 1.0
    while f(xk + t * dk) > fk + 1.0e-4 * t * slope
      t /= 1.5
    end
    return t
  end

function bfgs_quasi_newton_armijo(nlp ::  ADNLPModel, x0 :: Vector{T}) where {T<:Real}
    xk = x0
    sk = xk
    n = length(x0)
    fk = obj(nlp, xk)
    gk = grad(nlp, xk)
    yk = gk
    Hk = I(n)
    gnorm = gnorm0 = norm(gk)
    k = 0
    fk_inf = -1.0e20
    start_time = time()
    elapsed_time = time() - start_time
    max_time = 5.
    max_eval = 1000
    @printf "%2s %9s %9s\n" "k" "fk" "||∇f(x)||"
    @printf "%2d %9.2e %9.1e\n" k fk gnorm
    while gnorm > 1.0e-6 + 1.0e-6 * gnorm0 && k < 100 && fk > fk_inf && elapsed_time < max_time && neval_obj(nlp) < max_eval
      if k != 0
        ρk = 1.0/dot(yk, sk)
        Hk = (I-ρk*sk*yk') * Hk * (I-ρk*yk*sk') + ρk*sk*sk' 
      end
      dk = - Hk * gk  
      t = armijo(xk, dk, fk, gk, x -> obj(nlp, x))
      sk .= .-xk
      xk += t * dk
      sk .+= xk
      fk = obj(nlp, xk)
      yk .= .-gk
      gk = grad(nlp, xk)
      yk .+= gk
      gnorm = norm(gk)
      if k == 0
          Hk = dot(yk, sk) / dot(yk, yk) * I
      end
      k += 1
      elapsed_time = time() - start_time
      @printf "%2d %9.2e %9.1e %7.1e \n" k fk gnorm t
    end
    k > 100 && error("max iter")
    fk < fk_inf && error("non borné inf")
    elapsed_time > max_time && error("max time")
    neval_obj(nlp) > max_eval && error("max eval")
    return xk
  end
  
  bfgs_quasi_newton_armijo(nlp :: ADNLPModel) = bfgs_quasi_newton_armijo(nlp, nlp.meta.x0)

  #Test with Himmelblau function
fH(x) = (x[2]+x[1].^2-11).^2+(x[1]+x[2].^2-7).^2
gH(x) = [4*x[1]*(x[2]+x[1]^2-11)+2*(x[1]+x[2]^2-7);
    2*(x[2]+x[1]^2-11)+4*x[2]*(x[1]+x[2]^2-7)]
HH(x) = [ 4*(x[2]+3*x[1]^2-11)+2*x[1] 4*(x[1]+x[2]) ;
    4*(x[1]+x[2]) 2*x[2]+4*(x[1]+3*x[2]^2-7)] ; 
x0 = [10., 20.]

nlp3 = ADNLPModel(fH, x0)
bfgs_quasi_newton_armijo(nlp3)