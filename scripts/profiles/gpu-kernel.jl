using SparseArrays, CUDA
CUDA.allowscalar(false)

# function mul_Ax!(
#   y::CUDA.CuVector{T},
#   A_rowPtr,
#   A_colVal,
#   A_nzVal::CUDA.CuVector{T},
#   x::CUDA.CuVector{T},
#   n,
# ) where {T <: Real}
#   function kernel(y::CuDeviceVector{T}, A_rowPtr, A_colVal, A_nzVal, x, n)
#     row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     if row ≤ n
#       row_start = A_rowPtr[row]
#       row_end = A_rowPtr[row + 1] - 1
#       sumv = zero(T)
#       for element = row_start: row_end
#         sumv += A_nzVal[element] * x[A_colVal[element]]
#       end
#       y[row] = sumv
#     end
#     return nothing
#   end

#   threads = min(n, 256)
#   blocks = ceil(Int, n / threads)
#   @cuda name = "diagq" threads = threads blocks = blocks kernel(
#     y,
#     A_rowPtr,
#     A_colVal,
#     A_nzVal,
#     x,
#     n,
#   )
#   return y
# end

  # offset = warpsize() ÷ 2
  # while offset > 0
  #     val += shfl_down_sync(0xffffffff, val, offset)
  #     offset ÷= 2
  # end

function warp_reduce(val)
    # Loop unrolling for warpsize = 32
  val += shfl_down_sync(0xffffffff, val, 16, 32)
  val += shfl_down_sync(0xffffffff, val, 8, 32)
  val += shfl_down_sync(0xffffffff, val, 4, 32)
  val += shfl_down_sync(0xffffffff, val, 2, 32)
  val += shfl_down_sync(0xffffffff, val, 1, 32)
  return val
end

# function mul_Ax!(
#   y::CUDA.CuVector{T},
#   A_rowPtr,
#   A_colVal,
#   A_nzVal::CUDA.CuVector{T},
#   x::CUDA.CuVector{T},
#   n,
# ) where {T <: Real}

#   function kernel(y::CuDeviceVector{T}, A_rowPtr, A_colVal, A_nzVal, x, n)
#     thread_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1
#     warp_id = div(thread_id, 32)
#     lane = thread_id % 32
#     row = warp_id + 1
#     sumv = zero(T)
#     if row ≤ n
#       row_start = A_rowPtr[row]
#       row_end = A_rowPtr[row + 1] - 1
#       sumv = zero(T)
#       for element = row_start + lane: 32: row_end
#         sumv += A_nzVal[element] * x[A_colVal[element]]
#       end
#     end
#     sumv = warp_reduce(sumv)
#     if lane == 0 && row ≤ n
#       y[row] = sumv
#     end
#     return nothing
#   end

#   threads = min(n, 256)
#   blocks = ceil(Int, n / threads)
#   @cuda name = "mul2" threads = threads blocks = blocks kernel(
#     y,
#     A_rowPtr,
#     A_colVal,
#     A_nzVal,
#     x,
#     n,
#   )
#   return y
# end

function mul_Ax!(
  y::CUDA.CuVector{T},
  A_rowPtr,
  A_colVal,
  A_nzVal::CUDA.CuVector{T},
  x::CUDA.CuVector{T},
  n,
) where {T <: Real}

  function kernel(y::CuDeviceVector{T}, A_rowPtr, A_colVal, A_nzVal, x, n)
    thread_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1
    warp_id = div(thread_id, 32)
    lane = thread_id % 32
    row = warp_id + 1
    vals = @cuDynamicSharedMem(T, 32)
    if row ≤ n
      row_start = A_rowPtr[row]
      row_end = A_rowPtr[row + 1] - 1
      vals[threadIdx().x] = zero(T)
      for element = row_start + lane: 32: row_end
        vals[threadIdx().x] += A_nzVal[element] * x[A_colVal[element]]
      end
      (lane < 16) && (vals[threadIdx().x] += vals[threadIdx().x + 16])
      (lane < 8) && (vals[threadIdx().x] += vals[threadIdx().x + 8])
      (lane < 4) && (vals[threadIdx().x] += vals[threadIdx().x + 4])
      (lane < 2) && (vals[threadIdx().x] += vals[threadIdx().x + 2])
      (lane < 1) && (vals[threadIdx().x] += vals[threadIdx().x + 1])
      if lane == 0
        y[row] += vals[threadIdx().x] 
      end
    end
    return nothing
  end

  threads = min(n, 256)
  blocks = ceil(Int, n / threads)
  @cuda name = "mul2" threads = threads blocks = blocks kernel(
    y,
    A_rowPtr,
    A_colVal,
    A_nzVal,
    x,
    n,
  )
  return y
end

mul2!(y, A::CUDA.CUSPARSE.CuSparseMatrixCSR, x) = mul_Ax!(y, A.rowPtr, A.colVal, A.nzVal, x, length(x))

n = 10000
T = Float64
A = sprand(T, n, n, 0.2)
Agpu = CUDA.CUSPARSE.CuSparseMatrixCSR(A)
x, y = CUDA.rand(T, n), CUDA.rand(T, n)
mul2!(y, Agpu, x)
y - Agpu * x