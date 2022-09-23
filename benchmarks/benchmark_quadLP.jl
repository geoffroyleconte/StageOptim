using QuadraticModels, QPSReader
using Quadmath, SparseArrays, DoubleFloats
using RipQP
using CSV
using SolverBenchmark
# path_pb = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\quadLP\data\MPS"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"
path_pb = "/home/gelecd/datasets/quad_problems"
save_path = "/home/gelecd/code/docGL/benchmarks/ripqp_paper"

function createQuadraticModel_T(qpdata; T = Float128, name="qp_pb")
    return QuadraticModel(convert(Array{T}, qpdata.c), qpdata.qrows, qpdata.qcols,
            convert(Array{T}, qpdata.qvals),
            Arows=qpdata.arows, Acols=qpdata.acols,
            Avals=convert(Array{T}, qpdata.avals),
            lcon=convert(Array{T}, qpdata.lcon),
            ucon=convert(Array{T}, qpdata.ucon),
            lvar=convert(Array{T}, qpdata.lvar),
            uvar=convert(Array{T}, qpdata.uvar),
            c0=T(qpdata.c0), x0 = zeros(T, length(qpdata.c)), name=name)
end

T = Float128
Tlow = Float64

# compile
path_pb_lp = "/home/gelecd/.julia/artifacts/545f8c5577a056981a21caf3f53bd7b59cf67410/optrove-netlib-lp-f83996fca937"
qm1 = createQuadraticModel_T(readqps(string(path_pb_lp, "/AFIRO.SIF")), T=T)

ripqp_multik2(qm; T = T, Tlow = Tlow) = ripqp(qm, 
  mode = :multi,
  early_multi_stop = false,
  sp = K2KrylovParams{Tlow}( # solve in Float64
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 100,
    mem = 100,
    preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
    ρ_min=1.0e-15,
    δ_min = 1.0e-15,
    atol_min = 1.0e-16,
    rtol_min = 1.0e-16,
  ),
    sp2 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 5,
    mem = 5,
    preconditioner = LDL(T = T, pos = :R, warm_start = true),
    ρ_min=T(1.0e-15),
    δ_min = T(1.0e-15),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(),
  solve_method2=PC(),
  itol = InputTol(T, max_iter = 700, max_time = 2000.0, max_iter1 = 100, ϵ_pdd1 = T(1.0e1),
    ϵ_rc1 = T(1.0e-6), ϵ_rb1 = T(1.0e-6)),
  display = true,
)
stats = ripqp_multik2(qm1)

ripqp_multik3(qm; T = T, Tlow = Tlow) = ripqp(qm, 
mode = :multi,
early_multi_stop = false,
sp = K2KrylovParams{Tlow}( # solve in Float64
  uplo = :U,
  kmethod=:gmres,
  form_mat = true,
  equilibrate = false,
  itmax = 50,
  mem = 50,
  preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
  ρ_min=1.0e-15,
  δ_min = 1.0e-15,
  atol_min = 1.0e-16,
  rtol_min = 1.0e-16,
),
  sp2 = K2KrylovParams{T}( # solve in Float128
  uplo = :U,
  kmethod=:gmres,
  form_mat = true,
  equilibrate = false,
  itmax = 5,
  mem = 5,
  preconditioner = LDL(T = T, pos = :R, warm_start = true),
  ρ_min=T(1.0e-15),
  δ_min = T(1.0e-15),
  atol_min = T(1.0e-16),
  rtol_min = T(1.0e-16),
),
solve_method=IPF(),
solve_method2=PC(),
itol = InputTol(T, max_iter = 700, max_time = 10000.0, max_iter1 = 200, ϵ_pdd1 = T(1.0e1),
  ϵ_rc1 = T(1.0e-5), ϵ_rb1 = T(1.0e-5)),
display = true,
)
stats = ripqp_multik3(qm1)

ripqp_multik4(qm; T = T, Tlow = Tlow) = ripqp(qm, 
  mode = :multi,
  early_multi_stop = false,
  sp = K2KrylovParams{Tlow}( # solve in Float64
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 100,
    mem = 100,
    preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
    ρ_min=1.0e-15,
    δ_min = 1.0e-15,
    atol_min = 1.0e-16,
    rtol_min = 1.0e-16,
  ),
    sp2 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 5,
    mem = 5,
    preconditioner = LDL(T = T, pos = :R, warm_start = true),
    ρ_min=T(1.0e-15),
    δ_min = T(1.0e-15),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(),
  solve_method2=PC(),
  itol = InputTol(T, max_iter = 700, max_time = 2000.0, max_iter1 = 200, ϵ_pdd1 = T(1.0e1),
    ϵ_rc1 = T(1.0e-6), ϵ_rb1 = T(1.0e-6)),
  display = true,
)
stats = ripqp_multik4(qm1)

ripqp_multik5(qm; T = T, Tlow = Tlow) = ripqp(qm, 
  mode = :multi,
  early_multi_stop = false,
  sp = K2KrylovParams{Tlow}( # solve in Float64
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 100,
    mem = 10,
    preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
    ρ_min=1.0e-15,
    δ_min = 1.0e-15,
    atol_min = 1.0e-15,
    rtol_min = 1.0e-15,
  ),
    sp2 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 5,
    mem = 5,
    preconditioner = LDL(T = T, pos = :R, warm_start = true),
    ρ_min=T(1.0e-15),
    δ_min = T(1.0e-15),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(),
  itol = InputTol(T, max_iter = 700, max_time = 2000.0, max_iter1 = 100),
  display = true,
)
stats = ripqp_multik5(qm1)

ripqp_multik6(qm; T = T, Tlow = Tlow) = ripqp(qm, 
  mode = :multi,
  early_multi_stop = false,
  sp = K2KrylovParams{Tlow}( # solve in Float64
    uplo = :U,
    kmethod=:gmresir,
    form_mat = true,
    equilibrate = false,
    itmax = 100,
    mem = 10,
    preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
    ρ_min=1.0e-15,
    δ_min = 1.0e-15,
    atol_min = 1.0e-16,
    rtol_min = 1.0e-16,
  ),
    sp2 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 5,
    mem = 5,
    preconditioner = LDL(T = T, pos = :R, warm_start = true),
    ρ_min=T(1.0e-15),
    δ_min = T(1.0e-15),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(),
  solve_method2=PC(),
  itol = InputTol(T, max_iter = 700, max_time = 2000.0, max_iter1 = 100, ϵ_pdd1 = T(1.0e1),
    ϵ_rc1 = T(1.0e-6), ϵ_rb1 = T(1.0e-6)),
  display = true,
)
stats = ripqp_multik6(qm1)

ripqp_multik7(qm; T = T, Tlow = Tlow) = ripqp(qm, 
  mode = :multi,
  early_multi_stop = false,
  sp = K2KrylovParams{Tlow}( # solve in Float64
    uplo = :U,
    kmethod=:gmresir,
    form_mat = true,
    equilibrate = false,
    itmax = 100,
    mem = 10,
    preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
    ρ_min=1.0e-15,
    δ_min = 1.0e-15,
    atol_min = 1.0e-16,
    rtol_min = 1.0e-16,
  ),
    sp2 = K2KrylovParams{T}( # solve in Float128
    uplo = :U,
    kmethod=:gmres,
    form_mat = true,
    equilibrate = false,
    itmax = 5,
    mem = 5,
    preconditioner = LDL(T = T, pos = :R, warm_start = true),
    ρ_min=T(1.0e-15),
    δ_min = T(1.0e-15),
    atol_min = T(1.0e-16),
    rtol_min = T(1.0e-16),
  ),
  solve_method=IPF(),
  solve_method2=PC(),
  itol = InputTol(T, max_iter = 700, max_time = 2000.0, max_iter1 = 200, ϵ_pdd1 = T(1.0e1),
    ϵ_rc1 = T(1.0e-6), ϵ_rb1 = T(1.0e-6)),
  display = true,
)

ripqp_multi(qm; T = T, Tlow = Tlow) = ripqp(
  qm,
  mode = :multi,
  early_multi_stop = false,
  sp = K2LDLParams{Tlow}(),
  itol = InputTol(T, max_iter = 700, max_time = 10000.0, max_iter1 = 100),
)
stats = ripqp_multi(qm1)

ripqp_mono(qm; T = T) = ripqp(qm, itol = InputTol(T, max_iter = 700, max_time = 10000.0, max_iter1 = 100))
stats = ripqp_mono(qm1)

function optimize_ripqp(path_pb :: String, ripqp_func :: Function, T::DataType)
  problems = [
    createQuadraticModel_T(readqps(string(path_pb, "/TMA_ME_presolved.mps")), T = T, name = "TMA_ME"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAlift_presolved.mps")), T = T, name = "GlcAlift"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAerWT_presolved.mps")), T = T, name = "GlcAerWT"),
  ]

  return solve_problems(ripqp_func, problems)
end

function save_quad_problems(file_path :: String, ripqp_func :: Function; path_pb :: String = path_pb, T = T)
  lp_stats = optimize_ripqp(path_pb, ripqp_func, T)
  CSV.write(string(file_path, "_quad.csv"), lp_stats)
  return Nothing
end

function optimize_ripqp_nops(path_pb :: String, ripqp_func :: Function, T::DataType)
  problems = [
    createQuadraticModel_T(readqps(string(path_pb, "/TMA_ME.mps")), T = T, name = "TMA_ME"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAlift.mps")), T = T, name = "GlcAlift"),
    createQuadraticModel_T(readqps(string(path_pb, "/GlcAerWT.mps")), T = T, name = "GlcAerWT"),
  ]

  return solve_problems(ripqp_func, problems)
end

function save_quad_problems_nops(file_path :: String, ripqp_func :: Function; path_pb :: String = path_pb, T = T)
  lp_stats = optimize_ripqp_nops(path_pb, ripqp_func, T)
  CSV.write(string(file_path, "_nops_quad.csv"), lp_stats)
  return Nothing
end

# save_quad_problems(string(save_path, "/ripqp_multik2"), ripqp_multik2, T = T)
save_quad_problems(string(save_path, "/ripqp_multik3"), ripqp_multik3, T = T)
# save_quad_problems(string(save_path, "/ripqp_multik4"), ripqp_multik4, T = T)
# save_quad_problems(string(save_path, "/ripqp_multik5"), ripqp_multik5, T = T)
# save_quad_problems(string(save_path, "/ripqp_multik6"), ripqp_multik6, T = T)
# save_quad_problems(string(save_path, "/ripqp_multik7"), ripqp_multik7, T = T)
save_quad_problems(string(save_path, "/ripqp_multi1"), ripqp_multi, T = T)
save_quad_problems(string(save_path, "/ripqp_mono1"), ripqp_mono, T = T)
save_quad_problems_nops(string(save_path, "/ripqp_multik3"), ripqp_multik3, T = T)

# multik1
# ripqp_multik(qm; T = T, Tlow = Tlow) = ripqp(qm, 
#   mode = :multi,
#   early_multi_stop = false,
#   sp = K2KrylovParams{Tlow}( # solve in Float64
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = true,
#     itmax = 100,
#     mem = 100,
#     preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#     ρ_min=1.0e-8,
#     δ_min = 1.0e-8,
#     atol_min = 1.0e-15,
#     rtol_min = 1.0e-15,
#   ),
#   sp2 = K2KrylovParams{T}( # solve in Float128
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = true,
#     itmax = 20,
#     mem = 20,
#     preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#     ρ_min=T(1.0e-15),
#     δ_min = T(1.0e-15),
#     atol_min = T(1.0e-14),
#     rtol_min = T(1.0e-14),
#   ),
#     sp3 = K2KrylovParams{T}( # solve in Float128
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = false,
#     itmax = 5,
#     mem = 5,
#     preconditioner = LDL(T = T, pos = :R, warm_start = true),
#     ρ_min=T(1.0e-15),
#     δ_min = T(1.0e-15),
#     atol_min = T(1.0e-16),
#     rtol_min = T(1.0e-16),
#   ),
#   solve_method=IPF(),
#   ps=true,
#   itol = InputTol(
#     T,
#     max_iter = 700,
#     max_time = 4000.0,
#     max_iter1 = 100,
#     max_iter2 = 200,
#   ),
# )

# best so far
# ripqp_multik3(qm; T = T, Tlow = Tlow) = ripqp(qm, 
# mode = :multi,
# early_multi_stop = false,
# sp = K2KrylovParams{Tlow}( # solve in Float64
#   uplo = :U,
#   kmethod=:gmres,
#   form_mat = true,
#   equilibrate = false,
#   itmax = 50,
#   mem = 50,
#   preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#   ρ_min=1.0e-15,
#   δ_min = 1.0e-15,
#   atol_min = 1.0e-16,
#   rtol_min = 1.0e-16,
# ),
#   sp2 = K2KrylovParams{T}( # solve in Float128
#   uplo = :U,
#   kmethod=:gmres,
#   form_mat = true,
#   equilibrate = false,
#   itmax = 5,
#   mem = 5,
#   preconditioner = LDL(T = T, pos = :R, warm_start = true),
#   ρ_min=T(1.0e-15),
#   δ_min = T(1.0e-15),
#   atol_min = T(1.0e-16),
#   rtol_min = T(1.0e-16),
# ),
# solve_method=IPF(),
# solve_method2=PC(),
# itol = InputTol(T, max_iter = 700, max_time = 2000.0, max_iter1 = 200, ϵ_pdd1 = T(1.0e1),
#   ϵ_rc1 = T(1.0e-5), ϵ_rb1 = T(1.0e-5)),
# display = true,
# )