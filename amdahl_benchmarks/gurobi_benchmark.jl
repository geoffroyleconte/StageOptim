using Gurobi
using QPSReader
using QuadraticModels
using NLPModels
using SolverTools
using SolverBenchmark
using DataFrames
using LinearAlgebra
using JLD2
using SparseArrays

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end




function optimizeGurobi(QM)
    SM = SlackModel(QM)
    env = Gurobi.Env()

    # -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier,
    # 3=concurrent, 4=deterministic concurrent, 5=deterministic concurrent simplex.
    setparam!(env, "Method", 2)

    # set presolve to 0
    #setparam!(env, "Presolve", 0)

    # no scaling
    #setparam!(env, "ScaleFlag", 0)

    #setparam!(env, "Crossover", 0)

    Aeq = jac(SM, SM.meta.x0)
    beq = SM.meta.lcon
    f = grad(SM, zeros(length(SM.meta.x0)))
    H = hess(SM, zeros(length(SM.meta.x0)))
    H = Matrix(Symmetric(H, :L))
    n,m = size(Aeq)
    model = gurobi_model(env; f = f, H = H,
                        Aeq = Aeq, beq = beq,
                        lb = SM.meta.lvar, ub = SM.meta.uvar)
     # run optimization
    optimize(model)

    # y with dual: b'*y   s.t. A'*y <= c and y >= 0
    y = zeros(n)
    for i=1:n
        y[i] = Gurobi.get_dblattrelement(model, "Pi", i)
    end

    s = zeros(m) # s_l - s_u
    for i=1:m
        s[i] = Gurobi.get_dblattrelement(model, "RC", i)
    end

    # outputs
    optim_info = get_optiminfo(model)
    if optim_info.status == :optimal
        status = :acceptable
    elseif optim_info.status == :iteration_limit
        status = :max_iter
    elseif optim_info.status == :unbounded
        status = :unbounded
    else
        status = :unknown
    end

    x = get_solution(model)
    stats = GenericExecutionStats(status, SM, solution = x,
                                  objective = get_objval(model),
                                  iter = Gurobi.get_intattr(model,"BarIterCount"),
                                  primal_feas = norm(Aeq * x - beq),
                                  dual_feas = norm(Aeq' * y - H*x + s - f),
                                  solver_specific = Dict(:multipliers => y),
                                  elapsed_time = optim_info.runtime)
    return stats
end

function optimizeGurobi(qpdata::QPSData)
    return optimizeGurobi(createQuadraticModel(qpdata))
end

path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
path_pb_qp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/marosmeszaros"
# path_pb = "/user/eleves/gleconte2017/Z/Documents/TFE/marosmeszaros"
#path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
pb2 = string(path_pb_qp, "/DUAL1.SIF")
qpdata2 = readqps(pb2);
stats2 = optimizeGurobi(qpdata2)  # compile code

function optimize_netlib_gurobi(path_netlib)
    problems = []
    i_max = 1000
    i = 1
    for file_name in readdir(path_pb)
        if file_name[end-3:end] == ".SIF" && !(file_name in ["80BAU3B.SIF" ; "BORE3D.SIF";
                                                        "CAPRI.SIF"; "CZPROB.SIF";
                                                        "ETAMACRO.SIF"; "FINNIS.SIF";
                                                        "FORPLAN.SIF"; "GREENBEA.SIF";
                                                        "GREENBEB.SIF"; "MAROS.SIF";
                                                        "NESM.SIF"; "PEROLD.SIF";
                                                         "PILOT-JA.SIF"; "PILOT-WE.SIF";
                                                         "PILOT.SIF"; "PILOT4.SIF";
                                                         "PILOT87.SIF"; "PILOTNOV.SIF";
                                                          "RECIPELP.SIF"; "SHELL.SIF";
                                                         "SIERRA.SIF"; "STAIR.SIF";
                                                         "STANDATA.SIF"; "STANDGUB.SIF";
                                                        "STANDMPS.SIF"; "TUFF.SIF";
                                                        "VTP-BASE.SIF"; "DTOC3.SIF";
                                                         "HS35MOD.SIF";"QBORE3D.SIF";
                                                        "QCAPRI.SIF"; "QETAMACR.SIF";
                                                          "QFORPLAN.SIF"; "QPCSTAIR.SIF";
                                                        "QPCSTAIR.SIF"; "QPILOTNO.SIF";
                                                        "QRECIPE.SIF"; "QSHELL.SIF";
                                                        "QSIERRA.SIF"; "QSTAIR.SIF";
                                                        "QSTANDAT.SIF"; "UBH1.SIF";
                                                        "YAO.SIF"]) # problems with fixed variables


            println(file_name)
            pb_i = string(path_pb, "/", file_name)
            if file_name in ["BLEND.SIF"; "DFL001.SIF"; "FORPLAN.SIF"; "GFRD-PNC.SIF"; "SIERRA.SIF";
                        "EXDATA.SIF"; "QFORPLAN.SIF"; "QGFRDXPN.SIF"; "VALUES.SIF"]
                qpdata_i = readqps(pb_i, mpsformat=:fixed)
            else
                qpdata_i = readqps(pb_i)
            end
            push!(problems, createQuadraticModel(qpdata_i, name=file_name[1:end-4]))

            if i==i_max
                break
            end
            i += 1
        end
    end

    return solve_problems(optimizeGurobi, problems)
end



save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/StageOptim/amdahl_benchmarks/results"
# save_path = "/user/eleves/gleconte2017/Z/Documents/TFE/results"
#save_path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\results"

problems_stats_lp =  optimize_netlib_gurobi(path_pb_lp)

file_lp = jldopen(string(save_path, "/gurobi_test.jld2"), "w")
file_lp["stats"] = problems_stats_lp
close(file_lp)

# problems_stats_qp =  optimize_netlib_gurobi(path_pb_qp)
#
# file_qp = jldopen(string(save_path, "/mehrotra_qp11.jld2"), "w")
# file_qp["stats"] = problems_stats_qp
# close(file_qp)
