using LinearAlgebra, CUDA

# demo kernel
function gpu_add3!(y, x)
  index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  @inbounds y[index] += x[index]
  return
end

N = 10000
a, b = CUDA.rand(N), CUDA.rand(N)
numblocks = ceil(Int, N/256)

@cuda threads=256 blocks=numblocks gpu_add3!(a, b)

function bench_gpu3!(y, x)
  numblocks = ceil(Int, length(y)/256)
  CUDA.@sync begin
      @cuda threads=256 blocks=numblocks gpu_add3!(y, x)
  end
end

@btime bench_gpu3!($a, $b)