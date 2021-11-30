include("data_definition.jl")
include(raw"C:\Users\Geoffroy Leconte\.julia\dev\RipQP\src\RipQP.jl")
using Plots, FastTransforms

wk = [(k1 + 1/2) * ρ1 / (2 * m) for k1=0:n-1]
πℓ = [(k2 + 1/2) * ρ1 / (2 * m) for k2=0:n-1]
F = zeros(Complex{T}, n, n)
FFT = zeros(T, m+1, m+1)
FFTim = zeros(Complex{T}, m+1, m+1)
pFFT = plan_nufft1(wk, πℓ, 1.0e-10)
function mulopfft4!(res::AbstractVector{T}, f, α, β, F, FFT, FFTim, pFFT, X, Y, ξs, ηs, Δx, Δy) where {T}
  formF!(F, f, X, Y, n)
  FastTransforms.mul!(FFTim, pFFT, F)
  FFT .= real.(FFTim)
  display(FFT .* 4 .* Δx .* Δy)
  formfft!(res, FFT, ξs, ηs, m, 4 * α * Δx * Δy, β)
  return res
end

using QuadraticModels, LinearOperators, SparseArrays
opA4 = LinearOperator(T, ncon, nvar, false, false, 
                     (res, f, α, β) -> mulopfft4!(res, f, α, β, F, FFT, FFTim, pFFT, X, Y, ξs, ηs, Δx, Δy))