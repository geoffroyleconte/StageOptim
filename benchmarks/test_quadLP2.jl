using QuadraticModels, QPSReader
using Quadmath, SparseArrays, DoubleFloats
path_pb = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\quadLP\data\MPS"
# path_pb = "/home/mgi.polymtl.ca/geleco/quad_optim/problems/quadLP/data/MPS"
# path_pb = "/home/gelecd/datasets/quad_problems"

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

# qps1 = readqps(string(path_pb, "CYCLE.SIF"))
# qps1 = readqps(string(path_pb, "/TMA_ME_presolved.mps"))
# qps1 = readqps(string(path_pb, "/TMA_ME.mps"))
qps1 = readqps(string(path_pb, "/GlcAlift_presolved.mps"))
# qps1 = readqps(string(path_pb, "/GlcAerWT_presolved.mps"))

# using HSL
using RipQP
# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
T = Float128
qm1 = createQuadraticModel_T(qps1, T = T) # create QuadraticModel in Double64
Tlow = Float64
# stats1 = ripqp(qm1, mode = :multi, Timulti = Tlow,
#   sp = K2KrylovParams{Tlow}( # solve in double precision
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat  = true,
#     equilibrate = false,
#     itmax = 100,
#     mem = 100,
#     preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#     ρ_min=1.0e-30,
#     δ_min = 1.0e-30,
#     atol_min = 1.0e-16,
#     rtol_min = 1.0e-16,
#   ),
#   sp3 = K2KrylovParams{T}( # solve in quadruple precision
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat  = true,
#     equilibrate = false,
#     itmax = 100,
#     mem = 100,
#     preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#     ρ_min=T(1.0e-30),
#     δ_min = T(1.0e-30),
#     atol_min = T(1.0e-30),
#     rtol_min = T(1.0e-30),
#   ),
#   solve_method=IPF(),
#   ps=true,
#   itol = InputTol(
#     T,
#     ϵ_rb = T(1e-40), # very small to see what residuals can be reached
#     max_iter = 400,
#     max_time = 70000.0,
#     max_iter64 = 250,
#   ),
#   display = true,
# )

# println(stats1)


# T = Double64
qm1 = createQuadraticModel_T(qps1, T = T) # create QuadraticModel in Float128
Tlow = Float64
# stats2 = ripqp(qm1 , mode = :multi, sp = K2LDLParams{Float64}(), early_multi_stop = false,
#                 itol = InputTol(T, #ϵ_rb = T(1e-40), # very small to see what residuals can be reached
#                   # ϵ_rb64 = 1e-20, # very small to see what residuals can be reached
#                   max_iter = 700, max_time = 70000.0, max_iter1 = 100))

stats1 = ripqp(qm1, 
  mode = :multi,
  early_multi_stop = false,
  # sp = K2LDLParams(regul = :hybrid, ρ_min=1.0e-10, δ_min = 1.0e-10), # solve in Float64
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
    atol_min = 1.0e-16,
    rtol_min = 1.0e-16,
    Tir = T,
  ),
  # sp2 = K2KrylovParams{T}( # solve in Float128
  #   uplo = :U,
  #   kmethod=:gmres,
  #   form_mat = true,
  #   equilibrate = true,
  #   itmax = 20,
  #   mem = 20,
  #   preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
  #   ρ_min=T(1.0e-15),
  #   δ_min = T(1.0e-15),
  #   atol_min = T(1.0e-14),
  #   rtol_min = T(1.0e-14),
  # ),
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
    switch_solve_method = true,
  ),
  solve_method=IPF(),
  ps=true,
  itol = InputTol(
    T,
    # ϵ_rb = T(1e-40), # very small to see what residuals can be reached
    # ϵ_rb64 = 1e-20, # very small to see what residuals can be reached
    max_iter = 700,
    max_time = 70000.0,
    max_iter1 = 100,
    ϵ_pdd1 = T(1.0e1),
    ϵ_rc1 = T(1.0e-6),
    ϵ_rb1 = T(1.0e-6),
    # max_iter2 = 300,
  ),
  display = true,
  scaling = true,
)

# GlcAerWT
# stats1 = ripqp(qm1,
#          mode = :multi,
#          early_multi_stop = false,
#          # sp = K2LDLParams(regul = :hybrid, ρ_min=1.0e-10, δ_min = 1.0e-10), # solve in Float64
#          sp = K2KrylovParams{Tlow}( # solve in Float64
#            uplo = :U,
#            kmethod=:gmres,
#            form_mat = true,
#            equilibrate = false,
#            itmax = 100,
#            mem = 100,
#            preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#            ρ_min=1.0e-9,
#            δ_min = 1.0e-9,
#            atol_min = 1.0e-15,
#            rtol_min = 1.0e-15,
#          ),
#          sp2 = K2KrylovParams{T}( # solve in Float128
#            uplo = :U,
#            kmethod=:gmres,
#            form_mat = true,
#            equilibrate = false,
#            itmax = 20,
#            mem = 20,
#            preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#            ρ_min=T(1.0e-15),
#            δ_min = T(1.0e-15),
#            atol_min = T(1.0e-14),
#            rtol_min = T(1.0e-14),
#          ),
#            sp3 = K2KrylovParams{T}( # solve in Float128
#            uplo = :U,
#            kmethod=:gmres,
#            form_mat = true,
#            equilibrate = false,
#            itmax = 5,
#            mem = 5,
#            preconditioner = LDL(T = T, pos = :R, warm_start = true),
#            ρ_min=T(1.0e-15),
#            δ_min = T(1.0e-15),
#            atol_min = T(1.0e-16),
#            rtol_min = T(1.0e-16),
#          ),
#          solve_method=IPF(),
#          ps=true,
#          itol = InputTol(
#            T,
#            ϵ_rb = T(1e-40), # very small to see what residuals can be reached
#            # ϵ_rb64 = 1e-20, # very small to see what residuals can be reached
#            max_iter = 700,
#            max_time = 70000.0,
#            max_iter1 = 200,
#            max_iter2 = 300,
#          ),
#          display = true,
#          scaling = true,
#        )


# GlcAlift                  
# stats1 = ripqp(qm1, 
#   mode = :multi,
#   early_multi_stop = false,
#   # sp = K2LDLParams(regul = :hybrid, ρ_min=1.0e-10, δ_min = 1.0e-10), # solve in Float64
#   sp = K2KrylovParams{Tlow}( # solve in Float64
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = false,
#     itmax = 50,
#     mem = 50,
#     preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#     ρ_min=1.0e-10,
#     δ_min = 1.0e-10,
#     atol_min = 1.0e-15,
#     rtol_min = 1.0e-15,
#   ),
#   sp2 = K2KrylovParams{T}( # solve in Float128
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = false,
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
#     # ϵ_rb = T(1e-40), # very small to see what residuals can be reached
#     # ϵ_rb64 = 1e-20, # very small to see what residuals can be reached
#     max_iter = 700,
#     max_time = 70000.0,
#     max_iter1 = 100,
#     max_iter2 = 300,
#   ),
#   display = true,
#   scaling = true,
# )


# TMA_ME multi
# stats1 = ripqp(qm1, 
#   mode = :multi,
#   early_multi_stop = false,
#   # sp = K2LDLParams(regul = :hybrid, ρ_min=1.0e-10, δ_min = 1.0e-10), # solve in Float64
#   sp = K2KrylovParams{Tlow}( # solve in Float64
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = false,
#     itmax = 50,
#     mem = 50,
#     preconditioner = LDL(T = Tlow, pos = :R, warm_start = true),
#     ρ_min=1.0e-10,
#     δ_min = 1.0e-10,
#     atol_min = 1.0e-15,
#     rtol_min = 1.0e-15,
#   ),
#   sp2 = K2KrylovParams{T}( # solve in Float128
#     uplo = :U,
#     kmethod=:gmres,
#     form_mat = true,
#     equilibrate = false,
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
#     # ϵ_rb = T(1e-40), # very small to see what residuals can be reached
#     # ϵ_rb64 = 1e-20, # very small to see what residuals can be reached
#     max_iter = 300,
#     max_time = 70000.0,
#     max_iter1 = 100,
#     max_iter2 = 200,
#   ),
#   display = true,
#   scaling = true,
# )
