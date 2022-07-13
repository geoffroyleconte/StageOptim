include("data_definition.jl")
# include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using RipQP
using Plots

stor_K = fill!(Vector{T}(undef, length(DarkHole) * length(Pupil)), 4 * Δx * Δy)
c_K = 0
for (ξj1, ηj2) in DarkHole
  for (xk1, yk2) in Pupil
    c_K += 1
    stor_K[c_K] *= cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) 
  end
end

function mulopfft_opt3!(res::AbstractVector{T}, f, Pupil, DarkHole, stor_K) where {T}
  cDH = 0
  c_K = 0
  for (ξj1, ηj2) in DarkHole
    cDH += 1
    cP = 0
    fftj1j2 = zero(T)
    for (xk1, yk2) in Pupil
      c_K += 1
      cP += 1
      fftj1j2 += stor_K[c_K] * f[cP]
    end
    res[cDH] = fftj1j2
  end
end

stor_Kt = fill!(Vector{T}(undef, length(DarkHole) * length(Pupil)), 4 * Δx * Δy)
c_Kt = 0
for (xk1, yk2) in Pupil
  for (ξj1, ηj2) in DarkHole
    c_Kt += 1
    stor_Kt[c_Kt] *= cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) 
  end
end

function mulopfftt_opt3!(res::AbstractVector{T}, h, Pupil, DarkHole, stor_Kt) where {T}
  cP = 0
  c_Kt = 0
  for (xk1, yk2) in Pupil
    cDH = 0
    cP += 1
    fftk1k2 = zero(T)
    for (ξj1, ηj2) in DarkHole # ξj1 are sorted
      c_Kt += 1
      cDH += 1
      fftk1k2 += stor_Kt[c_Kt] * h[cDH]
    end
    res[cP] = fftk1k2
  end
end

using QuadraticModels, LinearOperators, SparseArrays
opA2 = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f) -> mulopfft_opt3!(res, f, Pupil, DarkHole, stor_K),
                     (res, f) -> mulopfftt_opt3!(res, f, Pupil, DarkHole, stor_Kt))
qm2 = QuadraticModel(c, LinearOperator(spzeros(T, nvar, nvar)), A = opA2, 
                    lcon = lcon, ucon = ucon, lvar = lvar, uvar = uvar)

stats1 = RipQP.ripqp(qm2, sp = RipQP.K1KrylovParams(uplo=:L, kmethod = :minres, mem = 1000, ρ_min = 1.0e-10, δ_min = 1.0e-10,
                        atol_min = 1.0e-10, rtol_min= 1.0e-10),
                        solve_method=RipQP.IPF(γ = 0.1), scaling = false, history=false, ps=false,
                        itol = RipQP.InputTol(ϵ_pdd = 1.0e-8, max_time=600 * 60.))

data = zeros(n, n)
c_out = 1
for indexes in idxPupil
  data[indexes.idx[1], indexes.idx[2]] = stats1.solution[c_out]
  c_out += 1
end
# using DelimitedFiles
# save_path = "/home/mgi.polymtl.ca/geleco/git_workspace/docGL/scripts/mask_astro/results"
# writedlm(string(save_path, "/lentil.csv"), data, ',')
heatmap(X, Y, data)
# using FFTW
# function mulopfft2!(res, f, α, β, Pupil, DarkHole, X, Y, ξs, ηs, Δx, Δy)
#   FFT2 = FFTW.r2r(f, FFTW.REDFT00)
#   FFTW.r2r!(FFT2, FFTW.REDFT00)
#   for (ξj1, ηj2) in DarkHole
#     cDH += 1
#     cP = 0
#     for (xk1, yk2) in Pupil
#       cP += 1
#       fftj1j2 = 4 * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * f[cP] *Δx * Δy
#     end
#     if β == 0
#       res[cDH] = α * fftj1j2
#     else
#       res[cDH] = α * fftj1j2 + β * res[cDH]
#     end
#   end
# end