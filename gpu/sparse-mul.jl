using LinearAlgebra, CUDA

@inline function warp_reduce(val)
    val += shfl_down_sync(0xffffffff, val, 16, 32)
    val += shfl_down_sync(0xffffffff, val, 8, 32)
    val += shfl_down_sync(0xffffffff, val, 4, 32)
    val += shfl_down_sync(0xffffffff, val, 2, 32)
    val += shfl_down_sync(0xffffffff, val, 1, 32)
    return val
end

function mul2!(y, A, x)

  function kernel(y, A_rowptr, A_colval, A_nzval, nrows, x)
    T = eltype(x)
    thread_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    CUDA.assume(warpsize() == 32)
    warp_size = 32
    warp_id = thread_id ÷ warp_size
    i = warp_id + 1
    lane = (thread_id - 1) % warp_size
    sum = zero(T)
    if i <= nrows
      for k = (A_rowptr[i] + lane): warp_size: (A_rowptr[i + 1] - 1)
        sum += A_nzval[k] * x[A_colval[k]]
      end
    end
    sum = warp_reduce(sum)
    if lane == 0 && i <= nrows
      y[i] = sum
    end
    return nothing
  end

  # config = launch_configuration(kernel.fun)
  # N = min(length(x), config.threads)
  threads = min(length(x), 256)
  # blocks = ceil(Int, length(A.colVal)/threads)
  blocks = ceil(Int, length(y)/threads)
  @cuda name="mul2" threads=threads blocks=blocks kernel(y, A.rowPtr, A.colVal, A.nzVal, size(A, 1), x)
  return y2
end