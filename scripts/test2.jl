using LinearAlgebra, Test

struct operator{T, F, Ft} 
    prod! :: F
    tprod :: Ft
end

function operator(A :: AbstractArray{T}) where T 
    prod! = (res, u, α, β) -> mul!(res, A, u, α, β)
    tprod! = (res, u, α, β) -> mul!(res, transpose(A), u, α, β)
    F = typeof(prod!)
    Ft = typeof(tprod!)
    return operator{T, F, Ft}(prod!, tprod!)
end

import LinearAlgebra.mul!
function mul!(res::Vector{T}, op::operator{T}, u::Vector{T}, α::T, β::T) where T 
    op.prod!(res, u, α, β)
end
function mul!(res::Vector{T}, op::operator{T}, u::Vector{T}) where T 
    op.prod!(res, u, one(T), zero(T))
end
function mul2!(res::Vector{T}, op::operator{T}, u::Vector{T}) where T 
    op.prod!(res, u, 2*one(T), 2*one(T))
end

function test_allocs()
    n = 10
    A = rand(n, n)
    A_op1 = operator(A)
    res1, u1, α1, β1 = rand(n), rand(n), 2.0, 2.0
    
    ### compile
    A_op1.prod!(res1, u1, β1, α1)
    mul!(res1, A_op1, u1, β1, α1)  
    mul!(res1, A, u1, α1, β1) 
    mul2!(res1, A_op1, u1) 
    ###

    ### tests allocs
    allocs1 = @allocated A_op1.prod!(res1, u1, α1, β1)
    @test allocs1 == 0
    allocs2 = @allocated mul!(res1, A_op1, u1, α1, β1)
    @test allocs2 == 0
    allocs3 = @allocated mul!(res1, A, u1, α1, β1) 
    @test allocs3 == 0
    allocs4 = @allocated mul!(res1, A_op1, u1) 
    @test allocs4 == 0
    allocs5 = @allocated mul2!(res1, A_op1, u1) 
    @test allocs5 == 0
end

test_allocs()