using QuadraticModels
using LinearAlgebra
using NLPModels
using QPSReader
using Printf

function createQuadraticModel(qpdata, name_pb)
    # probleme du point initial
    x0 = zeros(length(qpdata.lvar))
    for i=1:length(x0)
        if qpdata.lvar[i] == -Inf && qpdata.uvar[i] == Inf
            x0[i] = 1.
        elseif qpdata.lvar[i] == -Inf && qpdata.uvar[i] != Inf
            x0[i] = qpdata.uvar[i] - 1.
        elseif qpdata.lvar[i] != -Inf && qpdata.uvar[i] == Inf
            x0[i] = qpdata.lvar[i] + 1.
        else
            x0[i] = (qpdata.lvar[i] + qpdata.uvar[i]) / 2
        end
    end
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, x0=x0, name=name_pb)
end

function displayQuadraticModel(QM)
    #println("A = ", Matrix(jac(QM, QM.meta.x0)))
    #println("Q = ", Matrix(hess(QM, QM.meta.x0)))
    println("lvar = ", QM.meta.lvar)
    println("uvar = ", QM.meta.uvar)
    println("x0 = ", QM.meta.x0)
    #println("lcon = ", QM.meta.lcon)
    #println("ucon = ", QM.meta.ucon)
end
