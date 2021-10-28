using LinearAlgebra, SparseArrays, BenchmarkTools, Base.Threads, PrettyTables
# julia --threads 4

function custom_mul!(y::Vector{T}, A::SparseMatrixCSC{T}, x::Vector{T}) where T <: Number
  A.m == A.n || error("A is not a square matrix!")
  for i = 1 : A.n
    tmp = zero(T)
    for j = A.colptr[i] : (A.colptr[i+1] - 1)
      tmp += A.nzval[j] * x[A.rowval[j]]
    end
    y[i] = tmp
  end
  return y
end

function custom_inbounds_simd_mul!(y::Vector{T}, A::SparseMatrixCSC{T}, x::Vector{T}) where T <: Number
  A.m == A.n || error("A is not a square matrix!")
  for i = 1 : A.n
    tmp = zero(T)
    @inbounds @simd for j = A.colptr[i] : (A.colptr[i+1] - 1)
      tmp += A.nzval[j] * x[A.rowval[j]]
    end
    @inbounds y[i] = tmp
  end
  return y
end

function threaded_mul!(y::Vector{T}, A::SparseMatrixCSC{T}, x::Vector{T}) where T <: Number
  A.m == A.n || error("A is not a square matrix!")
  @threads for i = 1 : A.n
    tmp = zero(T)
    @inbounds for j = A.colptr[i] : (A.colptr[i+1] - 1)
      tmp += A.nzval[j] * x[A.rowval[j]]
    end
    @inbounds y[i] = tmp
  end
  return y
end

function threaded_simd_mul!(y::Vector{T}, A::SparseMatrixCSC{T}, x::Vector{T}) where T <: Number
  A.m == A.n || error("A is not a square matrix!")
  @threads for i = 1 : A.n
    tmp = zero(T)
    @inbounds @simd for j = A.colptr[i] : (A.colptr[i+1] - 1)
      tmp += A.nzval[j] * x[A.rowval[j]]
    end
    @inbounds y[i] = tmp
  end
  return y
end

function benchmark_mul!_T(T::DataType, n; convert_int = false)
  A = sprand(T, n, n, 0.2)
  if convert_int
    colptr, rowval = A.colptr, A.rowval
    colptr32, rowval32 = convert(Vector{Int32}, colptr), convert(Vector{Int32}, rowval)
    A = SparseMatrixCSC(n, n, colptr32, rowval32, A.nzval)
  end
  res = rand(T, n)
  b = rand(T, n)

  b1 = @benchmark mul!($res, $A, $b)
  b2 = @benchmark custom_mul!($res, $A, $b)
  b3 = @benchmark custom_inbounds_simd_mul!($res, $A, $b)
  b4 = @benchmark threaded_mul!($res, $A, $b)
  b5 = @benchmark threaded_simd_mul!($res, $A, $b)

  return [b1, b2, b3, b4, b5]
end

function benchmarks_mul!(n)
  row_names = [string(mul!),
               string(custom_mul!),
               string(custom_inbounds_simd_mul!),
               string(threaded_mul!),
               string(threaded_simd_mul!)]

  benchs64 = benchmark_mul!_T(Float64, n)
  benchs32 = benchmark_mul!_T(Float32, n)
  benchs32wInt = benchmark_mul!_T(Float32, n, convert_int = true)

  data = zeros(5, 3)
  data[:, 1] = [mean(bi.times) for bi in benchs64]
  data[:, 2] = [mean(bi.times) for bi in benchs32]
  data[:, 3] = [mean(bi.times) for bi in benchs32wInt]
  return data, row_names
end

data, row_names = benchmarks_mul!(10000)
t = pretty_table(data; 
  header = ["Float64", "Float32", "Float32wInt32"],
  row_names= row_names,
  title = "benchmarks single/double precision",
  formatters = ft_printf("%5.2e"),
  )