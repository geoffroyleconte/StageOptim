using LinearAlgebra


mutable struct Idx2D{T}
  idx::T
end
import Base.isless
isless(a::Idx2D, b::Idx2D) = a.idx[2] < b.idx[2]

T = Float64
n = 5 # domain discretization parameter
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
ξs = [i * ρ1 / m for i=0:m]
ηs = [i * ρ1 / m for i=0:m]
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