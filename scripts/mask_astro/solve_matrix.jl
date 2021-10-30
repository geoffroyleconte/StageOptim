include("data_definition.jl")
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using Plots

F = zeros(T, n, n)
FFT = zeros(T, m+1, m+1)
gj1 = zeros(T, n)

function compute_gj!(gj1, X, ξj1, F, Δx, n, transpose)
  for k2 = 1:n
    gj1k2 = zero(T)
    for k1 = 1:n
      Fk1k2 = !transpose ? F[k1, k2] : F[k2, k1]
      gj1k2 += 2 * cos(2 * π * X[k1] * ξj1) * F[k1, k2] * Δx
    end
    gj1[k2] = gj1k2
  end
end

function mulopfft!(res::AbstractVector{T}, FFT, f, F, α, β, gj1, X, Y, Δx, Δy, ξs, ηs, idxPupil, idxDarkHole, transpose) where {T}
  cDH = 0
  F .= zero(T)
  for k=1:length(idxPupil)
    idx2D = idxPupil[k]
    F[idx2D.idx[1], idx2D.idx[2]] = f[k]
  end
  for j1=1:m+1
    compute_gj!(gj1, X, ξs[j1], F, Δx, n, transpose)
    for j2=1:m+1
      fftj1j2 = zero(T)
      ηj2 = ηs[j2]
      for k2=1:n
        fftj1j2 += 2 * cos(2 * π * Y[k2] * ηj2) * gj1[k2] * Δy
      end
      FFT[j1, j2] = fftj1j2
    end
  end
  for k=1:length(idxDarkHole)
    idx2D = idxDarkHole[k]
    if β == zero(T)
      res[k] = α * FFT[idx2D.idx[1], idx2D.idx[2]]
    else
      res[k] = α * FFT[idx2D.idx[1], idx2D.idx[2]] + β * res[k]
    end
  end
end

using QuadraticModels, LinearOperators, SparseArrays
opA = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft!(res, FFT, f, F, α, β, gj1, X, Y, Δx, Δy, ξs, ηs, idxPupil, idxDarkHole, false),
                     (res, f, α, β) -> mulopfft!(res, FFT, f, F, α, β, gj1, X, Y, Δx, Δy, ξs, ηs, idxPupil, idxDarkHole, true))
qm = QuadraticModel(c, LinearOperator(spzeros(T, nvar, nvar)), A = opA, 
                    lcon = lcon, ucon = ucon, lvar = lvar, uvar = uvar)

stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
                     sp = RipQP.K1KrylovParams(uplo=:L, kmethod = :cg, preconditioner=:Identity),
                     solve_method=:IPF, scaling = false, history=false, presolve=false),
                     itol = RipQP.InputTol(ϵ_pdd = 1.0e-4, max_time=60.))

data = zeros(n, n)
c = 1
for indexes in idxPupil
  data[indexes.idx[1], indexes.idx[2]] = stats1.solution[c]
  c += 1
end
heatmap(X, Y, data)