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

function mulopfft!(res, f, α, β, Pupil, DarkHole)
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