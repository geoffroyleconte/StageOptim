using LinearAlgebra


mutable struct Idx2D{T}
  idx::T
end
import Base.isless
isless(a::Idx2D, b::Idx2D) = a.idx[2] < b.idx[2]

T = Float64
n = 500 # domain discretization parameter
ρ0 = T(4)
ρ1 = T(20)
sqrtϵ = T(1.0e-5)
Δx = T(1 / (2 * n))
Δy = Δx
X = [i * Δx for i=0.5:1:(n-0.5)]
Y = [i * Δy for i=0.5:1:(n-0.5)]
Pupil = [[x, y] for x in X  for y in Y if x^2 + y^2 < T(0.25)]

idxPupil = [Idx2D([i,j]) for i in 1:length(X)  for j in 1:length(Y) if X[i]^2 + Y[j]^2 < T(0.25)]
pv = sortperm(idxPupil)
# transposeidxPupil = idxPupil[perm][:]
nP = length(Pupil)
lvar = zeros(T, nP)
uvar = ones(T, nP)
# index fk1k2 = k1 + (k2-1) * n 

m = 35 # Fourier discretization parameter
Δξ = ρ1 / m
Δη = Δξ
ξs = [i * ρ1 / m for i=0:m]
ηs = [i * ρ1 / m for i=0:m]
check_DarkHole(ξ, η) = (ξ^2 + η^2 ≥ ρ0^2 && ξ^2 + η^2 ≤ ρ1^2 && η ≤ ξ)
DarkHole = [[ξ, η] for ξ in ξs for η in ηs if (ξ^2 + η^2 ≥ ρ0^2 && ξ^2 + η^2 ≤ ρ1^2 && η ≤ ξ)]
idxDarkHole = [Idx2D([i, j]) for i in 1:length(ξs) for j in 1:length(ηs) if (ξs[i]^2 + ηs[j]^2 ≥ ρ0^2 && ξs[i]^2 + ηs[j]^2 ≤ ρ1^2 && ηs[j] ≤ ξs[i])]
nDH = length(DarkHole)
# FFT = zeros(T, nDH, nP)
# cDH = 0
# for (ξj1, ηj2) in DarkHole
#   cDH += 1
#   cP = 0
#   for (xk1, yk2) in Pupil
#     fftj1j2 = 4 * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * Δx * Δy
#     cP += 1
#     FFT[cDH, cP] = fftj1j2
#   end
# end
lcon = fill(-sqrtϵ, nDH)
ucon = fill(sqrtϵ, nDH)
nvar, ncon = nP, nDH
c = ones(T, nvar) .* (-Δx * Δy)


function formFFT!(FFT, fft, ξs, ηs, m)
  FFT .= 0
  cFFT = 0
  for j2=1:m+1
    for j1=1:m+1
      if check_DarkHole(ξs[j1], ηs[j2])
        cFFT += 1
        FFT[j1, j2] = fft[cFFT]
      end
    end
  end
  return FFT
end
formFFT(fft, ξs, ηs, m) = formFFT!(zeros(m+1, m+1), fft, ξs, ηs, m)

function formF!(F, f, X, Y, n)
  F .= 0
  cF = 0
  for k2=1:n
    for k1=1:n
      if X[k1]^2 + Y[k2]^2 < T(0.25)
        cF += 1
        F[k1, k2] = f[cF]
      end
    end
  end
  return F
end
formF(f, X, Y, n) = formF!(zeros(n, n), f, X, Y, n)

function formf!(f, F, X, Y, n, α, β)
  f .= 0
  cF = 0
  for k2=1:n
    for k1=1:n
      if X[k1]^2 + Y[k2]^2 < T(0.25)
        cF += 1
        if β == 0
          f[cF] = α * F[k1, k2]
        else
          f[cF] = α * F[k1, k2] + β * f[cF]
        end
      end
    end
  end
  return f
end
formf(F, X, Y, n) = formf!(zeros(nP), F, X, Y, n, 1.0, 0.0)

function formfft!(fft, FFT, ξs, ηs, m, α, β)
  fft .= 0
  cFFT = 0
  for j2=1:m+1
    for j1=1:m+1
      if check_DarkHole(ξs[j1], ηs[j2])
        cFFT += 1
        if β == 0
          fft[cFFT] = α * FFT[j1, j2]
        else
          fft[cFFT] = α * FFT[j1, j2] + β * fft[cFFT]
        end
      end
    end
  end
  return fft
end
formfft(FFT, ξs, ηs, m) = formfft!(zeros(nDH), FFT, ξs, ηs, m, 1.0, 0.0)

# Pupil and DarkHole by column
idPupil = [[k1,k2] for k2 in 1:length(Y)  for k1 in 1:length(X) if X[k1]^2 + Y[k2]^2 < T(0.25)]
idDarkHole = [[j1, j2] for j2 in 1:length(ηs) for j1 in 1:length(ξs) if (ξs[j1]^2 + ηs[j2]^2 ≥ ρ0^2 && ξs[j1]^2 + ηs[j2]^2 ≤ ρ1^2 && ηs[j2] ≤ ξs[j1])]

# idx of the beginning of each column in idPupil and idDarkHole
colptr_DarkHole = ones(Int, m+1)
function get_colptr(idvec, n)
  nv = length(idvec)
  colptr = ones(Int, n+1)
  cv = 1
  col = 0
  for j=1:n
    k2 = idvec[cv][2]
    col = 0
    while k2 == j && cv < nv
      cv += 1
      k2 = idvec[cv][2]
      col += 1
    end
    if cv == nv && k2 == j
      col += 1
    end
    colptr[j+1] = col + colptr[j]
  end
  return colptr
end
colptr_Pupil = get_colptr(idPupil, n)
colptr_DarkHole = get_colptr(idDarkHole, m+1)
