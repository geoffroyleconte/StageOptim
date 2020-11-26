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

# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\quadLP\\data\\MPS"
path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"
QM = QuadraticModel(readqps(string(path_pb, "/TMA_ME.mps")))

env = CPLEX.Env()
CPXsetintparam(env, CPXPARAM_ScreenOutput, 1)   # Enable output (0=off)
CPXsetdblparam(env, CPXPARAM_TimeLimit, 3600)  # Time limit
# CPXsetintparam(env, CPXPARAM_Threads, 1) # Single thread
# CPXsetintparam(env, CPXPARAM_Barrier_Limits_Iteration, 0)
CPXsetintparam(env, CPXPARAM_SolutionType, 2)  # 2 = no crossover
CPXsetintparam(env, CPXPARAM_LPMethod, 4)  # 4 = Use barrier
CPXsetintparam(env, CPXPARAM_QPMethod, 4) # 4 = Use barrier, 0 = automatic
# Model m("theModel.mod");
# m.setExportFile("theModel.lp");
# m.stateConstraint();

status_p = Ref{Cint}()
lp = CPXcreateprob(env, status_p, "")
CPXnewcols(env, lp, QM.meta.nvar, QM.data.c, QM.meta.lvar, QM.meta.uvar, C_NULL, C_NULL)
CPXchgobjoffset(env, lp, QM.data.c0)
if QM.meta.nnzh > 0
    Hvals = zeros(eltype(QM.data.Hvals), length(QM.data.Hvals))
    for i=1:length(QM.data.Hvals)
        if QM.data.Hrows[i] == QM.data.Hcols[i]
            Hvals[i] = QM.data.Hvals[i] / 2
        else
            Hvals[i] = QM.data.Hvals[i]
        end
    end
    Q = sparse(QM.data.Hrows, QM.data.Hcols, QM.data.Hvals, QM.meta.nvar, QM.meta.nvar)
    diag_matrix = spdiagm(0 => diag(Q))
    Q = Q + Q' - diag_matrix
    qmatcnt = zeros(Int, QM.meta.nvar)
    for k = 1:QM.meta.nvar
      qmatcnt[k] = Q.colptr[k+1] - Q.colptr[k]
    end
    CPXcopyquad(env, lp, convert(Array{Cint,1}, Q.colptr[1:end-1].-1), convert(Array{Cint,1}, qmatcnt),
                convert(Array{Cint,1}, Q.rowval.-1), Q.nzval)
end

Acsrrowptr, Acsrcolval, Acsrnzval = sparse_csr(QM.data.Arows, QM.data.Acols,
                                               QM.data.Avals, QM.meta.ncon,
                                               QM.meta.nvar)

sense = fill(Cchar('A'), QM.meta.ncon) # lower, greater, range or equal. A is for the init
if length(QM.meta.jinf) > 0
    error("infeasible bounds in A")
end
p_low, p_upp, p_rng, p_fix = 1, 1, 1, 1
for j=1:QM.meta.ncon
   if length(QM.meta.jlow) > 0 && QM.meta.jlow[p_low] == j
       sense[j] = Cchar('G')
       if (p_low < length(QM.meta.jlow)) p_low += 1 end
   elseif length(QM.meta.jupp) > 0 && QM.meta.jupp[p_upp] == j
       sense[j] = Cchar('L')
       if (p_upp < length(QM.meta.jupp)) p_upp += 1 end
   elseif length(QM.meta.jrng) > 0 && QM.meta.jrng[p_rng] == j
       sense[j] = Cchar('R')
       if (p_rng < length(QM.meta.jrng)) p_rng += 1 end
   elseif length(QM.meta.jfix) > 0 && QM.meta.jfix[p_fix] == j
       sense[j] = Cchar('E')
       if (p_fix < length(QM.meta.jfix)) p_fix += 1 end
   else
       error("A error")
   end
end
rhs = zeros(QM.meta.ncon)
drange = zeros(QM.meta.ncon)
for j = 1:QM.meta.ncon
   if QM.meta.lcon[j] != -Inf && QM.meta.ucon[j] != Inf
       rhs[j] = QM.meta.ucon[j]
       drange[j] = QM.meta.ucon[j] - QM.meta.lcon[j]
   elseif QM.meta.lcon[j] != -Inf && QM.meta.ucon[j] == Inf
       rhs[j] = QM.meta.lcon[j]
   elseif QM.meta.lcon[j] == -Inf && QM.meta.ucon[j] != Inf
       rhs[j] = QM.meta.ucon[j]
   else
       rhs[j] = Inf
   end
end
CPXaddrows(env, lp, 0, QM.meta.ncon, length(Acsrcolval), rhs,
           sense, convert(Vector{Cint}, Acsrrowptr.- Cint(1)), convert(Vector{Cint}, Acsrcolval.- Cint(1)),
           Acsrnzval, C_NULL, C_NULL)

# save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\results\\cplex\\presolved_data"
save_path = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"
objoff = Ref{Cdouble}()
status = CPXpreslvwrite(env, lp, string(save_path, "/TMA_ME_presolved.pre"), objoff)

status_p = Ref{Cint}()
lp = CPXcreateprob(env, status_p, "")

status = CPXreadcopyprob(env, lp, string(save_path, "/TMA_ME_presolved.pre"), C_NULL)
status = CPXwriteprob(env, lp, string(save_path, "/TMA_ME_presolved.mps"), C_NULL)
