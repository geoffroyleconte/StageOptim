using BenchmarkTools, LinearAlgebra, SparseArrays
p = 8
time_Float32 = zeros(p)
time_Float64 = zeros(p)

function broadcast_test(z, x, y, α)
  z .= x .+ α .* y
end

for i = 1:p
  n = 10^i
  println("Dimension: ", n)
  T = Float32
  x = CUDA.rand(T, n)
  y = CUDA.rand(T, n)
  z = CUDA.rand(T, n)
  α = T(2.3)
  # @btime for k in 1:1000 dot($x, $y) end
  time_Float32[i] = @belapsed for k in 1:1000 broadcast_test($z, $x, $y, $α) end
end

for i = 1:p
  n = 10^i
  println("Dimension: ", n)
  T = Float64
  x = CUDA.rand(T, n)
  y = CUDA.rand(T, n)
  z = CUDA.rand(T, n)
  α = T(2.3)
  # @btime for k in 1:1000 dot($x, $y) end
  time_Float64[i] = @belapsed for k in 1:1000 broadcast_test($z, $x, $y, $α) end
end

using PrettyTables

dimension = [10^i for i=1:p]
time = hcat(dimension, time_Float32, time_Float64)
pretty_table(time ; header = ["Dimension", "Float32", "Float64"])

for i = 1:p
  val = min(time[i,2], time[i,3])
  time[i,2] /= val
  time[i,3] /= val
end
pretty_table(time ; header = ["Dimension", "Float32", "Flaot64"])