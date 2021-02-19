using QPSReader, SolverTools
using QuadraticModels
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
path_pb_QP = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"

mutable struct qps_presolve_data{T<:Real}
    activerows :: Vector{Bool}
    activecols :: Vector{Bool}
    rowcnt     :: Vector{Int} # number of nonzeros per row
    colcnt     :: Vector{Int}
    i_rm       :: Vector{Int}
    i_keep     :: Vector{Int}
    j_rm       :: Vector{Int}
    j_keep     :: Vector{Int}
    x          :: Vector{T}
end

function init_qpsp(qps)
    T = eltype(qps.c)
    qpsp = qps_presolve_data(fill(true, qps.ncon),
                             fill(true, qps.nvar),
                             zeros(Int, qps.ncon),
                             zeros(Int, qps.nvar),
                             Int[],
                             Int[],
                             Int[],
                             Int[],
                             zeros(T, qps.nvar)
                             )
    return qpsp
end

function empty_and_singleton_rows!(qps, qpsp)
    for arowi in qps.arows
        qpsp.rowcnt[arowi] += 1
    end
    for i=1:qps.ncon
        if qpsp.rowcnt[i] == 1 || qpsp.rowcnt[i] == 0
            qpsp.activerows[i] = false
        end
    end

    r_rm = 0 # nb of removed rows
    j = 0
    T= eltype(qps.c)
    for i=1:length(qps.arows)
        if !qpsp.activerows[qps.arows[i]]
            if qps.avals[i] > zero(T)
                qps.lvar[qps.acols[i]] = max(qps.lvar[qps.acols[i]], qps.lcon[qps.arows[i]] / qps.avals[i])
                qps.uvar[qps.acols[i]] = min(qps.uvar[qps.acols[i]], qps.ucon[qps.arows[i]] / qps.avals[i])
            elseif qps.avals[i] < zero(T)
                qps.lvar[qps.acols[i]] = max(qps.lvar[qps.acols[i]], qps.ucon[qps.arows[i]] / qps.avals[i])
                qps.uvar[qps.acols[i]] = min(qps.uvar[qps.acols[i]], qps.lcon[qps.arows[i]] / qps.avals[i])
            end
            @assert qps.lvar[qps.acols[i]] <= qps.uvar[qps.acols[i]]
        end
    end
end

function remove_rows!(qps, qpsp)
    qpsp.j_rm = findall(x->isequal(false, x), qpsp.activerows)
    qpsp.j_keep = zeros(Int, qps.ncon-length(qpsp.j_rm))
    j = 1
    for i=1:qps.ncon
        if qpsp.activerows[i] == true
            qpsp.j_keep[j] = i
            j += 1
        end
    end

    n_rm = length(qpsp.j_rm)
    if n_rm > 0
        j = 1
        n_max = length(qps.arows)
        for i=1:n_max
            r_rm = 0
            rm_row = false
            if qps.arows[i] >= qpsp.j_rm[1]
                r_rm = 1
                if qps.arows[i] > qpsp.j_rm[end]
                    r_rm = n_rm
                else
                    for k=1:n_rm
                        if qps.arows[i] == qpsp.j_rm[k]
                            rm_row = true
                            break
                        elseif qps.arows[i] < qpsp.j_rm[k]
                            r_rm = k-1
                            break
                        end
                    end
                end
            end
            qps.arows[j], qps.acols[j], qps.avals[j] = qps.arows[i]-r_rm, qps.acols[i], qps.avals[i]
            if !rm_row
                j += 1
            end
        end
        if j <= n_max
            resize!(qps.arows, j-1)
            resize!(qps.acols, j-1)
            resize!(qps.avals, j-1)
            sizehint!(qps.arows, j-1)
            sizehint!(qps.acols, j-1)
            sizehint!(qps.avals, j-1)
        end
        qps.lcon, qps.ucon = qps.lcon[qpsp.j_keep], qps.ucon[qpsp.j_keep]
        qps.ncon = length(qpsp.j_keep)
    end
end

function remove_ifix!(qps, qpsp) # detect ifix, change c, c0
    if length(qps.qrows) > 0
        T = eltype(qps.c)
        diag_Q = zeros(T, qps.nvar)
        for i=1:length(qps.qvals)
            if qps.qrows[i] == qps.qcols[i]
                diag_Q[qps.qrows[i]] = qps.qvals[i]
            end
        end
    end

    for i=1:qps.nvar
        if qps.lvar[i] == qps.uvar[i]
            qpsp.x[i] = qps.lvar[i]
            qpsp.activecols[i] = false
            if length(qps.qrows) > 0
                qps.c0 += qps.c[i] * qps.lvar[i] + diag_Q[i] * qps.lvar[i]^2 / 2
                for j=1:qps.nvar
                    qji = zero(T)
                    for k=1:length(qps.qvals)
                        if (qps.qrows[k] == i && qps.qcols[k] == j) || (qps.qrows[k] == j && qps.qcols[k] == i)
                            qji += qps.qvals[k]
                        end
                    end
                    qps.c[j] += qps.lvar[i] * qji
                end
            else
                qps.c0 += qps.c[i] * qps.lvar[i]
            end
        end
    end
end

function empty_and_singleton_cols!(qps, qpsp)
    if length(qps.qrows) > 0
        T = eltype(qps.c)
        for acoli in qps.acols
            qpsp.colcnt[acoli] += 1
        end
        n = length(qps.arows)
        for i=1:qps.nvar
            if qpsp.colcnt[i] == 0
                qpsp.activecols[i] = false
            elseif qpsp.colcnt[i] == 1 # slow, maybe worth to sort arows
                qpsp.activecols[i] = false
                idx = findfirst(qps.acols .== i)
                inf_constr, sup_constr = qps.lcon[qps.arows[idx]], qps.ucon[qps.arows[idx]]
                for k=1:n
                    if qps.acols[k] != qps.acols[idx] && qps.arows[k] == qps.arows[idx]
                        if qps.avals[k] > zero(T)
                            inf_constr -= qps.avals[k] * qps.uvar[qps.acols[k]]
                            sup_constr -= qps.avals[k] * qps.lvar[qps.acols[k]]
                        else
                            inf_constr -= qps.avals[k] * qps.lvar[qps.acols[k]]
                            sup_constr -= qps.avals[k] * qps.uvar[qps.acols[k]]
                        end
                    end
                end
                if qps.avals[idx] > 0
                    qps.lvar[idx] = max(qps.lvar[idx], qps.inf_constr / qps.avals[idx])
                    qps.uvar[idx] = min(qps.uvar[idx], qps.sup_constr / qps.avals[idx])
                elseif qps.avals[idx] < 0
                    qps.lvar[idx] = max(qps.lvar[idx], qps.sup_constr / qps.avals[idx])
                    qps.uvar[idx] = min(qps.uvar[idx], qps.inf_constr / qps.avals[idx])
                end
            end
            if qps.c[i] > 0
                qps.c0 += qps.c[i] * qps.lvar[i]
                qpsp.x[i] = qps.lvar[i]
            else
                qps.c0 += qps.c[i] * uvar
                qpsp.x[i] = qps.uvar[i]
            end
        end
    end
end

function remove_cols!(qps, qpsp)
    @assert issorted(qps.acols) && issorted(qps.qcols)
    # qpsp.i_keep = findall(x->isequal(true, x), qpsp.activecols)
    qpsp.i_rm = findall(x->isequal(false, x), qpsp.activecols)
    qpsp.i_keep = zeros(Int, qps.nvar-length(qpsp.i_rm))
    j = 1
    for i=1:qps.nvar
        if qpsp.activecols[i] == true
            qpsp.i_keep[j] = i
            j += 1
        end
    end
    n_rm = length(qpsp.i_rm)
    if n_rm > 0
        ## remove cols in A
        rm_idx = 1 # current index in i_rm
        n_max = length(qps.arows)  # 1st index to keep, last index to remove
        j = 1
        c_rm = 0 # nb of cols to remove
        max_idx = false
        for i in 1:n_max
            while rm_idx < n_rm && qps.acols[i] > qpsp.i_rm[rm_idx]
                rm_idx += 1
                c_rm += 1
            end
            if !max_idx && rm_idx == n_rm && qps.acols[i] > qpsp.i_rm[rm_idx]
                c_rm += 1
                max_idx = true
            end
            qps.arows[j], qps.acols[j], qps.avals[j] = qps.arows[i], qps.acols[i]-c_rm, qps.avals[i]
            if qps.acols[i] != qpsp.i_rm[rm_idx]
                j += 1
            else
                if qps.lvar[qps.acols[i]] == qps.uvar[qps.acols[i]]
                    qps.lcon[qps.arows[i]] -= qps.avals[i] * qps.lvar[qps.acols[i]]
                    qps.ucon[qps.arows[i]] -= qps.avals[i] * qps.lvar[qps.acols[i]]
                end
            end
        end
        if j <= n_max
            resize!(qps.arows, j-1)
            resize!(qps.acols, j-1)
            resize!(qps.avals, j-1)
            sizehint!(qps.arows, j-1)
            sizehint!(qps.acols, j-1)
            sizehint!(qps.avals, j-1)
        end
        ## remove cols in Q
        rm_idx = 1 # current index in i_rm
        n_max = length(qps.qrows)  # 1st index to keep, last index to remove
        j = 1
        c_rm = 0 # nb of cols to remove
        max_idx = false
        if n_max > 0
            for i in 1:n_max
                r_rm = 0
                rm_row = false
                if qps.qrows[i] >= qpsp.i_rm[1]
                    r_rm = 1
                    if qps.qrows[i] > qpsp.i_rm[end]
                        r_rm = n_rm
                    else
                        for k=1:n_rm
                            if qps.qrows[i] == qpsp.i_rm[k]
                                rm_row = true
                                break
                            elseif qps.qrows[i] < qpsp.i_rm[k]
                                r_rm = k-1
                                break
                            end
                        end
                    end
                end
                while rm_idx < n_rm && qps.qcols[i] > qpsp.i_rm[rm_idx]
                    rm_idx += 1
                    c_rm += 1
                end
                if !max_idx && rm_idx == n_rm && qps.qcols[i] > qpsp.i_rm[rm_idx]
                    c_rm += 1
                    max_idx = true
                end
                qps.qrows[j], qps.qcols[j], qps.qvals[j] = qps.qrows[i]-r_rm, qps.qcols[i]-c_rm, qps.qvals[i]
                if qps.qcols[i] != qpsp.i_rm[rm_idx] && !rm_row
                    j += 1
                end
            end
            if j <= n_max
                resize!(qps.qrows, j-1)
                resize!(qps.qcols, j-1)
                resize!(qps.qvals, j-1)
                sizehint!(qps.qrows, j-1)
                sizehint!(qps.qcols, j-1)
                sizehint!(qps.qvals, j-1)
            end
        end
        # qpsp.x_rm = qps.lvar[qpsp.i_rm]
        qps.c, qps.lvar, qps.uvar = qps.c[qpsp.i_keep], qps.lvar[qpsp.i_keep], qps.uvar[qpsp.i_keep]
        qps.nvar = length(qpsp.i_keep)
    end
end

function postsolve(qpsp, stats, QM)
    qpsp.x[qpsp.i_keep] = stats.solution
    stats_postsolve = GenericExecutionStats(stats.status, QM, solution = qpsp.x,
                                            objective = stats.objective,
                                            dual_feas = stats.dual_feas,
                                            primal_feas = stats.primal_feas,
                                            multipliers = stats.multipliers,
                                            multipliers_L = stats.multipliers_L,
                                            multipliers_U = stats.multipliers_U,
                                            iter = stats.iter,
                                            elapsed_time = stats.elapsed_time,
                                            solver_specific = stats.solver_specific)
    return stats_postsolve
end

function presolve!(qps)
    t = @timed begin
        qpsp = init_qpsp(qps)
        empty_and_singleton_rows!(qps, qpsp)
        remove_rows!(qps, qpsp)
        remove_ifix!(qps, qpsp)
        remove_cols!(qps, qpsp)
    end
    return qps, qpsp, t[2]
end

str = string(path_pb_QP, "\\QRECIPE.SIF")
qps = readqps(str)
qps, qpsp, t = presolve!(qps)

qm = QuadraticModel(qps)
stats1 = RipQP.ripqp(qm, mode=:mono, max_iter=400, regul=:classic, K=0, display=true)
println(stats1)

statsp1 = postsolve(qpsp, stats1, qm)
println(statsp1)


using SparseArrays
qptest = QPSData()
Q = sparse([1. 3. 0. 2.;
            0. 0. 0. 5.;
            0. 0. 4. 0.;
            0. 0. 0. 6.])
A = sparse([0. 0. 0. 1;
            0. 0. 0. 0;
            2. 3. 4. 4;
            0. 0. 0. 0;
            0. 0. 0. 0;
            0. 4. 0. 0;
            0. 2. 2. 0;
            0. 0. 0. 4;
            7. 8. 9. 0;
            0. 0. 0. 0;
            0. 0. 1. 0])
qptest.qrows, qptest.qcols, qptest.qvals = findnz(Q)
qptest.arows, qptest.acols, qptest.avals = findnz(A)
qptest.lvar = [-1.; -2.; 18; -3]
qptest.uvar = [-1.; -1.; 19.; -2]
qptest.c = [1. ;2.; 3.; 1.]
qptest.c0 = 1.
qptest.lcon= [0.; 0. ;1. ; 0.; 0.; -40.; 2.; 0.; -2.; 0.; 0.]
qptest.ucon= [0.; 0. ;5. ; 0.; 0.; 50.; 6.; 0.; 3.; 0.; 0.]
qptest.nvar, qptest.ncon = length(qptest.lvar), length(qptest.lcon)
qptest, qpsptest = presolve!(qptest)

display(Matrix(sparse(qptest.arows, qptest.acols, qptest.avals)))
# display(Matrix(sparse(qptest.qrows, qptest.qcols, qptest.qvals)))
