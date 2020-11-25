using CPLEX
using QuadraticModels, QPSReader
using SparseArrays, LinearAlgebra

function sparse_csr(I, J, V, m=maximum(I), n=maximum(J))
    csrrowptr = zeros(Int, m+1)
    # Compute the CSR form's row counts and store them shifted forward by one in csrrowptr
    coolen = length(I)
    min(length(J), length(V)) >= coolen || throw(ArgumentError("J and V need length >= length(I) = $coolen"))
    @inbounds for k in 1:coolen
        Ik = I[k]
        if 1 > Ik || m < Ik
            throw(ArgumentError("row indices I[k] must satisfy 1 <= I[k] <= m"))
        end
        csrrowptr[Ik+1] += 1
    end

    # Compute the CSR form's rowptrs and store them shifted forward by one in csrrowptr
    countsum = 1
    csrrowptr[1] = 1
    @inbounds for i in 2:(m+1)
        overwritten = csrrowptr[i]
        csrrowptr[i] = countsum
        countsum += overwritten
    end

    # Counting-sort the column and nonzero values from J and V into csrcolval and csrnzval
    # Tracking write positions in csrrowptr corrects the row pointers
    csrcolval = zeros(Int, length(I))
    csrnzval = zeros(length(I))
    @inbounds for k in 1:coolen
        Ik, Jk = I[k], J[k]
        if 1 > Jk || n < Jk
            throw(ArgumentError("column indices J[k] must satisfy 1 <= J[k] <= n"))
        end
        csrk = csrrowptr[Ik+1]
        csrrowptr[Ik+1] = csrk + 1
        csrcolval[csrk] = Jk
        csrnzval[csrk] = V[k]
    end
    csrrowptr = csrrowptr[1:end-1]
    return csrrowptr, csrcolval, csrnzval
end


env = CPLEX.Env()
CPXsetintparam(env, CPXPARAM_ScreenOutput, 1)   # Enable output (0=off)
CPXsetdblparam(env, CPXPARAM_TimeLimit, 3600)  # Time limit
CPXsetintparam(env, CPXPARAM_Threads, 1) # Single thread
for (k, v) in kwargs
    if k==:presolve
        CPXsetintparam(env, CPXPARAM_Preprocessing_Presolve, v) # 0 = no presolve
    elseif k==:scaling
        CPXsetintparam(env, CPXPARAM_Read_Scale, -1) # -1 = no scaling
    elseif k==:crossover
        CPXsetintparam(env, CPXPARAM_SolutionType, 2)  # 2 = no crossover
    end
end
CPXsetintparam(env, CPXPARAM_LPMethod, method)  # 4 = Use barrier
CPXsetintparam(env, CPXPARAM_QPMethod, method) # 4 = Use barrier, 0 = automatic
