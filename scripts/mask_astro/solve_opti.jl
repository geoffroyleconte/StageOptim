include("data_definition.jl")
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using Plots

function mulopfft!(res::AbstractVector{T}, f, α, β, Pupil, DarkHole, Δx, Δy) where {T}
  cDH = 0
  for (ξj1, ηj2) in DarkHole
    cDH += 1
    cP = 0
    fftj1j2 = zero(T)
    for (xk1, yk2) in Pupil
      cP += 1
      fftj1j2 += 4 * α * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * f[cP] * Δx * Δy
    end
    if β == 0
      res[cDH] = fftj1j2
    else
      res[cDH] = fftj1j2 + β * res[cDH]
    end
  end
end

gj1 = zeros(T, n)

function mulopfftt!(res::AbstractVector{T}, h, α, β, Pupil, DarkHole, Δx, Δy) where {T}
  cP = 0
  for (xk1, yk2) in Pupil
    cDH = 0
    cP += 1
    fftk1k2 = zero(T)
    for (ξj1, ηj2) in DarkHole # ξj1 are sorted
      cDH += 1
      fftk1k2 += 4 * α * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * h[cDH] * Δx * Δy
    end
    if β == 0
      res[cP] = fftk1k2
    else
      res[cP] = fftk1k2 + β * res[cP]
    end
  end
end

using QuadraticModels, LinearOperators, SparseArrays
opA = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft!(res, f, α, β, Pupil, DarkHole, Δx, Δy),
                     (res, f, α, β) -> mulopfftt!(res, f, α, β, Pupil, DarkHole, Δx, Δy))
qm = QuadraticModel(c, LinearOperator(spzeros(T, nvar, nvar)), A = opA, 
                    lcon = lcon, ucon = ucon, lvar = lvar, uvar = uvar)

stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                        sp = RipQP.K1KrylovParams(uplo=:L, kmethod = :cg, preconditioner=:Identity),
                        solve_method=:IPF, scaling = false, history=false, presolve=false),
                        itol = RipQP.InputTol(ϵ_pdd = 1.0e-4, max_time=5 * 60.))

data = zeros(n, n)
c = 1
for indexes in idxPupil
  data[indexes.idx[1], indexes.idx[2]] = stats1.solution[c]
  c += 1
end
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