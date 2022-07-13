# max tr(C X) subject to X is positive semidefinite
                      #  tr(A_i X) = b_i for i = 1, ...,m
using LinearAlgebra
C = [1 0 0 0;
     0 2 0 0;
     0 0 3 0;
     0 0 0 4]

A1 = [1 0 0 0;
      0 1 0 0;
      0 0 0 0;
      0 0 0 0]

A2 = [0 0 0 0;
      0 1 0 0;
      0 0 5 2;
      0 0 2 6]

b = [10.0; 20.0]

using ProxSDP, Convex
n = 4
X = Semidefinite(n)
problem = maximize(dot(C, X), dot(A1, X) == b[1], dot(A2, X) == b[2])

# Solve optimization problem with ProxSDP
solve!(problem, ProxSDP.Optimizer(log_verbose=true, tol_gap=1e-4, tol_feasibility=1e-4))

# Get the objective value
problem.optval # -30

# Retrieve solution
evaluate(X)