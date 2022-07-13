include("data_definition.jl")
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using Plots

F = zeros(T, n, n)
FFT = zeros(T, m+1, m+1)
gj1 = zeros(T, n)
K = zeros(T, m+1, n)
for j1=1:m+1
  for k2=1:n
    K[j1, k2] = 2 * cos(2 * π * X[k2] * ξs[j1]) * Δx
  end
end
Ki = zeros(T, n, m+1)
for k1=1:n
  for j2=1:m+1
    Ki[k1, j2] = 2 * cos(2 * π * X[k1] * ξs[j2]) * Δξ
  end
end

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

function mulopfft!(res::AbstractVector{T}, FFT, f, F, α, β, X, Y, Δx, Δy, ξs, ηs) where {T}
  cDH = 0
  F .= zero(T)
  cF = 0
  for k2=1:n
    for k1=1:n
      if X[k1]^2 + Y[k2]^2 < T(0.25)
        cF += 1
        F[k1, k2] = f[cF]
      end
    end
  end
  for j1=1:m+1
    ξj1 = ξs[j1]
    for j2=1:m+1
      ηj2 = ηs[j2]
      fftj1j2 = zero(T)
      if check_DarkHole(ξj1, ηj2)
        for k1=1:n
          xk1 = X[k1]
          for k2=1:n
            yk2 = Y[k2]
            fftj1j2 += 4 * cos(2 * π * xk1 * ξj1) * cos(2 * π * yk2 * ηj2) * F[k1, k2] * Δy * Δx
          end
        end
      end
      FFT[j1, j2] = fftj1j2
    end
  end
  cFFT = 0
  for j2=1:m+1
    for j1=1:m+1
      if check_DarkHole(ξs[j1], ηs[j2])
        cFFT += 1
        if β == 0
          res[cFFT] = α * FFT[j1, j2]
        else
          res[cFFT] = α * FFT[j1, j2] + β * res[cFFT]
        end
      end
    end
  end
end


function mulopFFT!(FFT::AbstractMatrix{T}, F, X, Y, Δx, Δy, ξs, ηs) where {T}
  for j1=1:m+1
    ξj1 = ξs[j1]
    for j2=1:m+1
      ηj2 = ηs[j2]
      fftj1j2 = zero(T)
      if check_DarkHole(ξj1, ηj2)
        for k1=1:n
          xk1 = X[k1]
          for k2=1:n
            yk2 = Y[k2]
            fftj1j2 += 4 * cos(2 * π * xk1 * ξj1) * cos(2 * π * yk2 * ηj2) * F[k1, k2] * Δy * Δx
          end
        end
      end
      FFT[j1, j2] = fftj1j2
    end
  end
  return FFT
end

function mulopfftt!(res::AbstractVector{T}, FFT, fft, F, α, β, X, Y, Δξ, Δη, ξs, ηs) where {T}
  cDH = 0
  FFT .= zero(T)
  F .= zero(T)
  cFFT = 0
  for j2=1:m+1
    for j1=1:m+1
      if check_DarkHole(ξs[j1], ηs[j2])
        cFFT += 1
        FFT[j1, j2] = fft[cFFT]
      end
    end
  end
  for k1=1:n
    xk1 = X[k1]
    for k2=1:n
      yk2 = Y[k2]
      fk1k2 = zero(T)
      if xk1^2 + yk2^2 < T(0.25)
        for j1=1:m+1
          ξj1 = ξs[j1]
          for j2=1:m+1
            ηj2 = ηs[j2]
            fk1k2 += 4 * cos(2 * π * xk1 * ξj1) * cos(2 * π * yk2 * ηj2) * FFT[j1, j2] * Δξ * Δη
          end
        end
      end
      F[k1, k2] = fk1k2
    end
  end
  cF = 0
  for k2=1:n
    for k1=1:n
      if X[k1]^2 + Y[k2]^2 < T(0.25)
        cF += 1
        if β == 0
          res[cF] = α * F[k1, k2]
        else
          res[cF] = α * F[k1, k2] + β * res[cF]
        end
      end
    end
  end
end

tmpM, FFT, F = zeros(T, n, m+1), zeros(T, m+1, m+1), zeros(T, n, n)
function mulopFFTmat!(fft::AbstractVector{T}, K, f, α, β, F, FFT, tmpM, X, Y, ξs, ηs, n, m) where {T}
  formF!(F, f, X, Y, n)
  tmpM .= 0
  mul!(tmpM, F, K')
  mul!(FFT, K, tmpM)
  formfft!(fft, FFT, ξs, ηs, m, α, β)
end

function mulopFFTTmat!(f::AbstractVector{T}, K, fft, α, β, F, FFT, tmpM, X, Y, ξs, ηs, n, m) where {T}
  formFFT!(FFT, fft, ξs, ηs, m)
  tmpM .= 0
  mul!(tmpM, K', FFT)
  mul!(F, tmpM, K)
  formf!(f, F, X, Y, n, α, β)
end

using QuadraticModels, LinearOperators, SparseArrays
opAmat = LinearOperator(T, ncon, nvar, false, false, 
                        (fft, f, α, β) -> mulopFFTmat!(fft, K, f, α, β, F, FFT, tmpM, X, Y, ξs, ηs, n, m),
                        (f, fft, α, β) -> mulopFFTTmat!(f, K, fft, α, β, F, FFT, tmpM, X, Y, ξs, ηs, n, m))

#prod!
a = rand(nP)
Fa = formF(a, X, Y, n)
opA = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft!(res, FFT, f, F, α, β, X, Y, Δx, Δy, ξs, ηs),
                     (res, fft, α, β) -> mulopfftt!(res, FFT, fft, F, α, β, X, Y, Δξ, Δη, ξs, ηs))
afft = opA * a
aFFT = formFFT(afft, ξs, ηs, m)
aFFT3 = mulopFFT!(zeros(m+1, m+1), Fa, X, Y, Δx, Δy, ξs, ηs)
aFFT2 = K * Fa * K'
for i=1:m+1
  for j=1:m+1
    if !check_DarkHole(ξs[i], ηs[j])
      aFFT2[i,j] = 0.
    end
  end
end

#tprod!
function mulopFFTT!(F::AbstractMatrix{T}, FFT, X, Y, Δξ, Δη, ξs, ηs) where {T}
  for k1=1:n
    xk1 = X[k1]
    for k2=1:n
      yk2 = Y[k2]
      fk1k2 = zero(T)
      if xk1^2 + yk2^2 < T(0.25)
        for j1=1:m+1
          ξj1 = ξs[j1]
          for j2=1:m+1
            ηj2 = ηs[j2]
            fk1k2 += 4 * cos(2 * π * xk1 * ξj1) * cos(2 * π * yk2 * ηj2) * FFT[j1, j2] * Δξ * Δη
          end
        end
      end
      F[k1, k2] = fk1k2
    end
  end
  return F
end
Fb = mulopFFTT!(zeros(n, n), aFFT2, X, Y, Δξ, Δη, ξs, ηs)


qm = QuadraticModel(c, LinearOperator(spzeros(T, nvar, nvar)), A = opAmat, 
                    lcon = lcon, ucon = ucon, lvar = lvar, uvar = uvar)

stats1 = RipQP.ripqp(qm, sp = RipQP.K1KrylovParams(uplo=:L, kmethod = :minres),
                     solve_method=RipQP.IPF(), scaling = false, history=false, ps=false,
                     itol = RipQP.InputTol(ϵ_pdd = 1.0e-4, max_time=60., max_iter=30))

datap = zeros(n, n)
c = 1
for indexes in idxPupil
  datap[indexes.idx[1], indexes.idx[2]] = stats1.solution[c]
  c += 1
end
heatmap(X, Y, datap)

