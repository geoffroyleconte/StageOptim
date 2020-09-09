using Tulip
using QPSReader, SparseArrays, QuadraticModels, LinearAlgebra
using SolverTools, SolverBenchmark
using JLD2

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end

function tulip_solve(QM)
    model = Tulip.Model{Float64}()
    Tulip.set_parameter(model, "BarrierIterationsLimit", 200)
    Tulip.set_parameter(model, "TimeLimit", 1200.)
    #Tulip.set_parameter(model, "Presolve", 0)
    Tulip.set_parameter(model, "BarrierTolerancePFeas", 1.0e-6)
    Tulip.set_parameter(model, "BarrierToleranceDFeas", 1.0e-6)
    A = sparse(QM.data.Arows, QM.data.Acols, QM.data.Avals, QM.meta.ncon, QM.meta.nvar)
    varnames = [string("X", i) for i=1:length(QM.data.c)]
    connames = [string("X", i) for i=1:length(QM.meta.lcon)]
    Tulip.load_problem!(model.pbdata,
        QM.meta.name,
        QM.meta.minimize, QM.data.c, QM.data.c0,
        A,
        QM.meta.lcon, QM.meta.ucon,
        QM.meta.lvar, QM.meta.uvar,
        connames, varnames
    )
    start_time = time()
    Tulip.optimize!(model)
    elapsed_time = time() - start_time
    # rbNorm =  norm(model.solution.Ax .- model.solver.b, Inf)
    y = model.solution.y_lower .- model.solution.y_upper
    s = model.solution.s_lower .- model.solution.s_upper
    # rcNorm = norm(A' * y .+ s .- dat.c, Inf)
    rbNorm = model.solver.res.rp_nrm
    rcNorm = model.solver.res.rd_nrm
    if model.solver.solver_status == Tulip.Trm_Optimal
        status = :acceptable
    elseif model.solver.solver_status == Trm_TimeLimit
        status = :max_time
    elseif model.solver.solver_status == Trm_IterationLimit
        status = :max_iter
    elseif model.solver.solver_status == Trm_PrimalInfeasible ||
                model.solver.solver_status == Trm_DualInfeasible ||
                model.solver.solver_status == Trm_PrimalDualInfeasible
        status = :infeasible
     elseif model.solver.solver_status == Trm_NumericalProblem ||
                model.solver.solver_status == Trm_MemoryLimit
        status = :exception
    else
        status = :unknown
    end

    stats = GenericExecutionStats(status, QM, solution = model.solution.x,
                                  objective = model.solution.z_primal,
                                  dual_feas = rcNorm,
                                  primal_feas = rbNorm,
                                  multipliers = y,
                                  multipliers_L =  model.solution.s_lower,
                                  multipliers_U = model.solution.s_upper,
                                  iter = model.solver.niter,
                                  elapsed_time=elapsed_time)
    return stats
end

# path_pb_lp = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
# dat = readqps(string(path_pb_lp, "\\AFIRO.SIF"))
path_pb_lp = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/netlib"
dat = readqps(string(path_pb_lp, "/AFIRO.SIF"))
QM = QuadraticModel(dat)
stats1 = tulip_solve(QM)


function optimize_tulip(path_pb)
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
    return solve_problems(tulip_solve, problems)
end

save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/StageOptim/amdahl_benchmarks/results"
problems_stats_lp =  optimize_tulip(path_pb_lp)

file_lp = jldopen(string(save_path, "/tulip_lp2.jld2"), "w")
file_lp["stats"] = problems_stats_lp
close(file_lp)
