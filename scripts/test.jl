using QPSReader, QuadraticModels, SolverCore, SolverBenchmark, SparseArrays, TimerOutputs, LinearAlgebra
using QDLDL
# using RipQP
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib_ps"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
# qm = QuadraticModel(readqps(string(path_pb, "\\BANDM_PS.mps")))
qm = QuadraticModel(readqps(string(path_pb, "\\AFIRO.SIF"), mpsformat=:fixed))
# stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(refinement = :none, kc=0,mode=:mono, scaling=true, 
#                      sp = RipQP.K2_5hybridParams(preconditioner = :ActiveCHybridLDL), solve_method=:PC),
#                      itol = RipQP.InputTol(max_iter=100, ϵ_rb32 = 1e-6) )#,
TimerOutputs.enable_debug_timings(RipQP)
reset_timer!(RipQP.to)
stats1 = RipQP.ripqp(qm, 
                    # sp = RipQP.K2LDLParams(safety_dist_bnd=true, fact_alg = RipQP.LDLFact(regul=:classic)),
                     # w = RipQP.SystemWrite(write=false, name=string(save_path, "/bug_minres"),kfirst=1, kgap=10),
                    #  sp = RipQP.K1CholParams(),
                    # sp = RipQP.K2KrylovParams(uplo = :L, kmethod = :minres, rhs_scale=false, #δ0 = 0.,
                    # form_mat = true,
                    #   equilibrate = false,
                    #   preconditioner = RipQP.LDL(fact_alg=RipQP.LDLFact(), T = Float32, warm_start = true, pos=:R),
                    #   ρ_min=1.0e-8, δ_min = 1.0e-8,
                    #   mem = 10,
                    #   itmax = 10,
                    #   atol_min = 1.0e-8, rtol_min = 1.0e-8,
                    #   # k3_resid = true,
                    #   # cb_only = true,
                    #   ),   
                      # sp2 = RipQP.K2KrylovParams(uplo = :U, kmethod = :gmres, rhs_scale=true, #δ0 = 0.,
                      # form_mat = true,
                      #   # equilibrate = false,
                      #   preconditioner = RipQP.LDL(fact_alg=RipQP.LDLFact(), T = Float64, pos = :R, warm_start = true),
                      #   ρ_min=1.0e-8, δ_min = 1.0e-8,
                      #   mem = 5,
                      #   itmax = 5,
                      #   atol0 = 1.0e-2, rtol0=1.0e-2,
                      #   atol_min = 1.0e-10, rtol_min = 1.0e-10,
                      #   # k3_resid = true,
                      #   # cb_only = true,
                      #   ), 
                    # sp2 = RipQP.K2LDLParams(fact_alg = RipQP.LDLFact(regul = :classic), bypass_bound_dist_safety = true),
                    #  sp = RipQP.K3SKrylovParams(uplo = :U, kmethod=:minres, rhs_scale=true, #δ0 = 0.,
                    #         preconditioner = RipQP.BlockDiagK3S(),
                    #         ρ_min=1.0e-8, δ_min = 1.0e-8,
                    #         mem = 100,
                    #         atol_min = 1.0e-6, rtol_min = 1.0e-6,
                    #         ), 
                    solve_method=RipQP.PC(), 
                    scaling = true, history=false, ps=true, mode=:mono, kc=0,
                    #  solve_method2=RipQP.PC(),
                     perturb = false,
                     # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                     itol = RipQP.InputTol(max_iter=400, max_time=100.0, max_iter1 = 30,
                      #  ϵ_rc=1.0e-6, ϵ_rb=1.0e-6, ϵ_pdd=1.0e-8,
                     ),
                     display = true,
                   )
println(stats1)
println(maximum(-(Symmetric(qm.data.H, :L) * stats1.solution) + (qm.data.A' * stats1.multipliers) + 
    stats1.multipliers_L - stats1.multipliers_U - qm.data.c))
# TimerOutputs.complement!(RipQP.to)
# show(RipQP.to, sortby = :firstexec)
println(sum(stats1.solver_specific[:KresNormH]))

print(aze+z)

using DelimitedFiles, MatrixMarket
K = MatrixMarket.mmread(string(save_path, "\\bug_minaresK_iter1.mtx"))
rhs_aff = readdlm(string(save_path, "\\bug_minaresrhs_iter1_aff.rhs"), Float64)[:]
rhs_cc =  readdlm(string(save_path, "\\bug_minaresrhs_iter1_cc.rhs"), Float64)[:] 

function riptest(qm)
    stats = RipQP.ripqp(qm, display = false, sp = RipQP.K2LDLParams(), solve_method=RipQP.PC(), scaling = true, 
    ps = true,
   itol = RipQP.InputTol(max_iter=150, max_time=10.0),
  #  ϵ_rc=1.0e-1, ϵ_rb=1.0e-1, ϵ_pdd=1.0e0,
   )
  #  println(maximum(-(Symmetric(qm.data.H, :L) * stats.solution) + (qm.data.A' * stats.multipliers) + 
  #   stats.multipliers_L - stats.multipliers_U - qm.data.c))
   return stats
end

function createQuadraticModel(qpdata; name="qp_pb")
    return QuadraticModel(qpdata.c, qpdata.qrows, qpdata.qcols, qpdata.qvals,
            Arows=qpdata.arows, Acols=qpdata.acols, Avals=qpdata.avals,
            lcon=qpdata.lcon, ucon=qpdata.ucon, lvar=qpdata.lvar, uvar=qpdata.uvar,
            c0=qpdata.c0, name=name)
end

problems = []
i_max = 14
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
        pb_i = string(path_pb, "\\", file_name)
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
problems_stats = solve_problems(riptest, problems,
                                colstats=[:name, :status, :elapsed_time, :objective, :dual_feas, :primal_feas, :iter]);


using DataFrames
using SolverBenchmark
using SolverTools
using JLD2
res_path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\StageOptim\\amdahl_benchmarks\\results"
fgur3 = jldopen(string(res_path, "\\gurobi_scaling_lp3.jld2"), "r")  # no display
gurobi3 = fgur3["stats"];
close(fgur3)
fcplex1 = jldopen(string(res_path, "\\cplex_scaling_lp3.jld2"), "r")
cplex1 = fcplex1["stats"];
close(fcplex1)
fxpress1 = jldopen(string(res_path, "\\xpress_scaling_lp2.jld2"), "r")
xpress1 = fxpress1["stats"];
close(fxpress1)
file_rip_mono1 = jldopen(string(res_path, "\\G-2021-03_qp_mono.jld2"), "r")
ripqp_mono1 = file_rip_mono1["stats"];
close(file_rip_mono1)
file_rip_mp4 = jldopen(string(res_path, "\\ripqp_qp_mp4.jld2"), "r")
ripqp_mp4 = file_rip_mp4["stats"];
close(file_rip_mp4)
# dynamic regul2 multi precision default params
file_rip_mp_d1 = jldopen(string(res_path, "\\G-2021-03_qp_mono_c.jld2"), "r")
ripqp_mp_d1 = file_rip_mp_d1["stats"];
close(file_rip_mp_d1)
# stats = Dict(:xpress=>xpress1,
#              :ripqp=>ripqp_mono1,
#              :cplex=>cplex1,
#              :gurobi=>gurobi3
#              )
stats = Dict(:ripqp=>ripqp_mono1,
             :ripqp_centrality_corr=>ripqp_mp_d1)

             
using Plots
perf = performance_profile(stats, df->df.iter)
plot!(perf, legend=:bottomright)
title!("Performance profile (Netlib problems)")
savefig(raw"C:\Users\Geoffroy Leconte\Documents\cours\TFE\code\results\performance_profile\profileG-03_net_iter_d.pdf")

# exemple création qm
# probleme1

# using QuadraticModels
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using LinearAlgebra, SparseArrays, SparseMatricesCOO
Q = [6. 2. 1.
     2. 5. 2.
     1. 2. 4.]
c = [-8.; -3; -3]
A = [1. 0. 1.
     0. 2. 1.]
b = [0.; 3]
l = [0.;0;0]
u = [Inf; Inf; Inf]
QM = QuadraticModel(c, SparseMatrixCOO(tril(Q)), A=SparseMatrixCOO(A), lcon=b, ucon=b, lvar=l, uvar=u, c0=0., name="QM1");
# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\QuadraticModels.jl\src\QuadraticModels.jl")
# QM = QuadraticModels.QuadraticModel(c, Q, A=A, lcon=[-3; -4], ucon=[-2.; Inf], lvar=l, uvar=u, c0=0., name="QM1");
stats1 = RipQP.ripqp(QM,
                sp = RipQP.K2KrylovParams(
                      form_mat = true,
                       uplo=:L, 
                       preconditioner = RipQP.Equilibration(), 
                       kmethod = :minres, ρ_min = 1e0 * sqrt(eps()), δ_min = 1e0 * sqrt(eps())),
                solve_method=RipQP.IPF(), scaling = false, history=false, ps=false, 
                itol = RipQP.InputTol(max_time = 20.0, max_iter = 50))
# stats1 = RipQP.ripqp(QM)
println(stats1)


#mLCP
using QPSReader, QuadraticModels, SolverTools, SolverBenchmark
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_marosmeszaros"
qm = QuadraticModel(readqps(string(path_pb, "\\KB2.SIF"), mpsformat=:fixed))
# probleme1
Q = [6. 2 1
     2 5 2
     1 2 4]
c = [-8.; -3; -3]
c0 = 0.
A = [1. 0 1
    0 2 1]
b = [0.; 3]
lvar = [0. ,0 ,0]
uvar = [Inf, Inf, Inf]
mLCP = RipQP.mLCPModel(c, -b, sparse(Q), sparse(A'),
                        sparse(A), sparse([0. 0; 0 0]), lvar, uvar)
stats1 = RipQP.ripmLCP(mLCP, display=true)
display(stats1)

function custom_mul!(res, M, x)
  res .= 0
  @inbounds for j ∈ axes(M, 2)
    @simd for k ∈ nzrange(M, j)
      i = M.rowval[k]
      res[i] += M.nzval[k] * x[j]
    end
  end
end

using DelimitedFiles, MatrixMarket, LinearAlgebra
K = MatrixMarket.mmread( string(save_path, "\\CVXQP1_MK_iter1.mtx"))
K = K + K' - Diagonal(K)
using AMD, Metis
p1 = amd(K)
K[p1,p1]
p2, ip2 = Metis.permutation(K)
K[p2,p2]
p3 = Metis.partition(K, 3)
K[p3, p3]

# Float16
function QuadraticModel_T(qm, T)
  return QuadraticModel(
    Vector{T}(qm.data.c),
    spzeros(T, qm.meta.nvar, qm.meta.nvar),
    A = T.(qm.data.A),
    lcon = Vector{T}(qm.meta.lcon),
    ucon = Vector{T}(qm.meta.ucon),
    lvar = Vector{T}(qm.meta.lvar),
    uvar = Vector{T}(qm.meta.uvar),
    x0 = zeros(T, qm.meta.nvar),
  )
end
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
# using RipQP
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_netlib_ps"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\problemes_marosmeszaros"
save_path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\systems"
# path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\doctorat\\code\\datasets\\lptestset"
# qm = QuadraticModel(readqps(string(path_pb, "\\BANDM_PS.mps")))
qm = QuadraticModel(readqps(string(path_pb, "\\AFIRO.SIF"), mpsformat=:fixed))
qm16 = QuadraticModel_T(qm, Float16)
# stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(refinement = :none, kc=0,mode=:mono, scaling=true, 
#                      sp = RipQP.K2_5hybridParams(preconditioner = :ActiveCHybridLDL), solve_method=:PC),
#                      itol = RipQP.InputTol(max_iter=100, ϵ_rb32 = 1e-6) )#,
stats1 = RipQP.ripqp(qm16, iconf = RipQP.InputConfig(
                        # w = RipQP.SystemWrite(write=false, name=string(save_path, "/bug_minres"),kfirst=1, kgap=10),
                        sp = RipQP.K2LDLParams(),
                        # sp = RipQP.K2KrylovParams(uplo = :L, kmethod=:gmres,# δ_min = 0.,
                        #         # ρ0=0., δ0 = 0.,
                        #         ), 
                        solve_method=:IPF, scaling = true, history=true, presolve=false,
                        # w = RipQP.SystemWrite(write=true, kfirst=1, name = string(save_path, "\\CVXQP1_M"), kgap=1000)), 
                        ),
                     itol = RipQP.InputTol(Float16, max_iter = 20),
                     display = true,
                     )
println(stats1)

using LDLFactorizations, SparseArrays, LinearAlgebra, BenchmarkTools
n = 1000
A = sprand(Float64, n, n, 0.1)
A = triu(A * A' + I)
A32 = convert(SparseMatrixCSC{Float32, Int}, A)
LDL = ldl(Symmetric(A, :U))
LDL32 = ldl(Symmetric(A32, :U))
res = rand(Float64, n);
res32 = rand(Float32, n);

@benchmark ldiv!($LDL, $res)
@benchmark ldiv!($LDL32, $res32)