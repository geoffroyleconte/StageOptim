include("data_definition.jl")
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using Plots

Xline = zeros(Int, n) # number of Y for each X
for i=1:n
  cx = 0
  for j=1:n
    if X[i]^2 + Y[j]^2 < T(0.25)
      cx += 1
    end
  end
  Xline[i] = cx
end
csXline = [zero(Int); cumsum(Xline)] 

ξsline = zeros(Int, m + 1)
for i=1:m+1
  cξ = 0
  for j=1:m+1
    if check_DarkHole(ξs[i], ηs[j])
      cξ += 1
    end
  end
  ξsline[i] = cξ
end
csξsline = [zero(Int); cumsum(ξsline)] 

# gⱼ₁ₖ₂ = 2 Σₖ₁ cos(2πxₖ₁ξⱼ₁) fₖ₁ₖ₂ Δx
function compute_gj1!(gj1, X, ξj1, f, Δx, n, Xline, csXline)
  cF = 0
  gj1 .= 0
  for k2 = 1:n
    for k1 = 1:n
      if X[k1]^2 + X[k2]^2 < T(0.25) 
        cF += 1
        gj1[k2] += 2 * cos(2 * π * X[k1] * ξj1) * f[cF] * Δx
      end
    end
  end
end

# gₖ₁ⱼ₂ = 2 Σₖ₂ cos(2πyₖ₂ηⱼ₂) fₖ₁ₖ₂ Δy
gj2 = zeros(T, n)
function compute_gj2!(gj2, X, Y, ηj2, f, Δy, n)
  cF = 0
  gj2 .= 0
  for k2 = 1:n
    for k1 = 1:n
      if X[k1]^2 + Y[k2]^2 < T(0.25) 
        cF += 1
        gj2[k1] += 2 * cos(2 * π * Y[k2] * ηj2) * f[cF] * Δy
      end
    end
  end
end

gj1 = zeros(T, n)
gk2 = zeros(T, m+1)
# fftⱼ₁ⱼ₂ = 2 Σₖ₂ cos(2π yₖ₂ ηⱼ₂) gⱼ₁ₖ₂ Δy
function mulopfft2!(res::AbstractVector{T}, f, α, β, X, Y, ξs, ηs, gj2, Δx, Δy) where {T}
  cDH = 0
  for j2=1:m+1
    ηj2 = ηs[j2]
    compute_gj2!(gj2, X, Y, ηj2, f, Δy, n)
    for j1=1:m+1
      ξj1 = ξs[j1]
      if check_DarkHole(ξj1, ηj2)
        cDH += 1
        fftj1j2 = zero(T)
        for k1=1:n
          fftj1j2 += 2 * α * cos(2*π*X[k1]*ξj1) * gj2[k1] * Δx
        end
        if β == 0
          res[cDH] = fftj1j2
        else
          res[cDH] = fftj1j2 + β * res[cDH]
        end
      end
    end
  end
end

# gₖ₁ⱼ₂ = 2 Σⱼ₁ cos(2πxₖ₁ξⱼ₁) hⱼ₁ⱼ₂ Δx
function compute_gk2!(gk2, ξs, ηs, yk2, fft, Δη, m)
  cFFT = 0
  gk2 .= 0
  for j2 = 1:m+1
    ηj2 = ηs[j2]
    for j1 = 1:m+1
      ξj1 = ξs[j1]
      if check_DarkHole(ξj1, ηj2)
        cFFT += 1
        gk2[j1] += 2 * cos(2 * π * yk2 * ηj2) * fft[cFFT] * Δη
      end
    end
  end
end

function mulopfftt2!(res::AbstractVector{T}, h, α, β, Pupil, DarkHole, Δx, Δy) where {T}
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

function mulopfftt3!(res::AbstractVector{T}, fft, α, β, X, Y, ξs, ηs, gk2, Δξ, Δη) where {T}
  cP = 0
  for k2=1:n
    yk2 = Y[k2]
    compute_gk2!(gk2, ξs, ηs, yk2, fft, Δη, m)
    for k1=1:n
      xk1 = X[k1]
      if xk1^2 + yk2^2 < T(0.25)
        cP += 1
        fk1k2 = zero(T)
        for j1=1:m+1
          fk1k2 += 2 * α * cos(2*π*xk1*ξs[j1]) * gk2[j1] * Δξ
        end
        if β == 0
          res[cP] = fk1k2
        else
          res[cP] = fk1k2 + β * res[cP]
        end
      end
    end
  end
end

using QuadraticModels, LinearOperators, SparseArrays
opA2 = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft2!(res, f, α, β, X, Y, ξs, ηs, gj2, Δx, Δy),
                     (res, f, α, β) -> mulopfftt2!(res, f, α, β, Pupil, DarkHole, Δx, Δy))

opA3 = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft2!(res, f, α, β, X, Y, ξs, ηs, gj2, Δx, Δy),
                     (res, fft, α, β) -> mulopfftt3!(res, fft, α, β, X, Y, ξs, ηs, gk2, Δx, Δy))

# qm = QuadraticModel(c, LinearOperator(spzeros(T, nvar, nvar)), A = opA3, 
#                     lcon = lcon, ucon = ucon, lvar = lvar, uvar = uvar)

# stats1 = RipQP.ripqp(qm, iconf = RipQP.InputConfig(
#                         sp = RipQP.K1KrylovParams(uplo=:L, kmethod = :cg, preconditioner=:Identity),
#                         solve_method=:IPF, scaling = false, history=false, presolve=false),
#                         itol = RipQP.InputTol(ϵ_pdd = 1.0e-4, max_time=10 * 60.))

# data_sol = zeros(n, n)
# cnt = 1
# for indexes in idxPupil
#   data_sol[indexes.idx[1], indexes.idx[2]] = stats1.solution[cnt]
#   cnt += 1
# end
# heatmap(X, Y, data_sol)
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