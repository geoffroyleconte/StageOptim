using LinearAlgebra
T = Float64
n = 10 # domain discretization parameter
ρ0 = T(4)
ρ1 = T(20)
ϵ = T(1.0e-3)
Δx = T(1 / (2 * n))
Δy = Δx
X = [i * Δx for i=0.5:1:(n-0.5)]
Y = [i * Δy for i=0.5:1:(n-0.5)]
Pupil = [(x, y) for x in X  for y in Y if x^2 + y^2 < T(0.25)]
idxPupil = [(i,j) for i in 1:length(X)  for j in 1:length(Y) if X[i]^2 + Y[j]^2 < T(0.25)]
nP = length(Pupil)
lvar = zeros(T, nP)
uvar = ones(T, nP)
# index fk1k2 = k1 + (k2-1) * n 

m = 35 # Fourier discretization parameter
ξs = [i * ρ1 / m for i=0:m]
ηs = [i * ρ1 / m for i=0:m]
DarkHole = [(ξ, η) for ξ in ξs for η in ηs if (ξ^2 + η^2 ≥ ρ0^2 && ξ^2 + η^2 ≤ ρ1^2 && η ≤ ξ)]
nDH = length(DarkHole)
FFT = zeros(T, nDH, nP)
cDH = 0
for (ξj1, ηj2) in DarkHole
  cDH += 1
  cP = 0
  for (xk1, yk2) in Pupil
    fftj1j2 = 4 * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * Δx * Δy
    cP += 1
    FFT[cDH, cP] = fftj1j2
  end
end
lcon = fill(-ϵ, nDH)
ucon = fill(ϵ, nDH)
nvar, ncon = nP, nDH
c = ones(T, nvar)

function mulopfft!(res::AbstractVector{T}, f, α, β, Pupil, DarkHole) where {T}
  cDH = 0
  for (ξj1, ηj2) in DarkHole
    cDH += 1
    cP = 0
    fftj1j2 = zero(T)
    for (xk1, yk2) in Pupil
      cP += 1
      fftj1j2 += 4 * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * f[cP] *Δx * Δy
    end
    if β == 0
      res[cDH] = α * fftj1j2
    else
      res[cDH] = α * fftj1j2 + β * res[cDH]
    end
  end
end

function mulopfftt!(res::AbstractVector{T}, g, α, β, Pupil, DarkHole) where {T}
  cP = 0
  for (xk1, yk2) in Pupil
    cDH += 0
    cP += 1
    fftk1k2 = zero(T)
    for (ξj1, ηj2) in DarkHole
      cDH += 1
      fftk1k2 += 4 * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * g[cDH] *Δx * Δy
    end
    if β == 0
      res[cDH] = α * fftk1k2
    else
      res[cDH] = α * fftk1k2 + β * res[cDH]
    end
  end
end

using QuadraticModels, LinearOperators, SparseArrays
opA = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft!(res, f, α, β, Pupil, DarkHole),
                     (res, g, α, β) -> mulopfftt!(res, g, α, β, Pupil, DarkHole))
qm = QuadraticModel(c, spzeros(T, nvar, nvar), A = opA, 
                    lcon = lcon, ucon = ucon, lvar = lvar, uvar = uvar)


using FFTW
function mulopfft2!(res, f, α, β, Pupil, DarkHole, X, Y, ξs, ηs, Δx, Δy)
  FFT2 = FFTW.r2r(f, FFTW.REDFT00)
  FFTW.r2r!(FFT2, FFTW.REDFT00)
  for (ξj1, ηj2) in DarkHole
    cDH += 1
    cP = 0
    for (xk1, yk2) in Pupil
      cP += 1
      fftj1j2 = 4 * cos(2*π*xk1*ξj1) * cos(2*π*yk2*ηj2) * f[cP] *Δx * Δy
    end
    if β == 0
      res[cDH] = α * fftj1j2
    else
      res[cDH] = α * fftj1j2 + β * res[cDH]
    end
  end
end