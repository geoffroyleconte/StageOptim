using LinearAlgebra, LDLFactorizations, SparseArrays

function create_K1_2(AT, LH, D, δ, id, qp)
    T = eltype(D)
    nnzK = length(AT.rowval) + length(LH.rowval) + id.nvar # A, LH, D
    K_rowval = Vector{Int}(undef, nnzK) 
    K_nzval = Vector{T}(undef, nnzK) 
    if qp
        K_colptr = Vector{Int}(undef, id.ncon+2*id.nvar+1)
        fill_K1_2!(K_colptr, K_rowval, K_nzval, AT.colptr, AT.rowval, AT.nzval, LH.colptr, LH.rowval, LH.nzval, 
                   D, δ, id.nvar, id.ncon)

        return SparseMatrixCSC(id.nvar, id.ncon+2*id.nvar, K_colptr, K_rowval, K_nzval)
    else
        K_colptr = Vector{Int}(undef, id.ncon+id.nvar+1)
        fill_K1_2lp!(K_colptr, K_rowval, K_nzval, AT.colptr, AT.rowval, AT.nzval, D, δ, id.nvar, id.ncon)

        return SparseMatrixCSC(id.nvar, id.ncon+id.nvar, K_colptr, K_rowval, K_nzval)
    end  
end

function fill_K1_2!(K_colptr, K_rowval, K_nzval, AT_colptr, AT_rowval, AT_nzval, LH_colptr, LH_rowval, LH_nzval, D, δ, nvar, ncon)

    # fill first bloc
    K_colptr[1:ncon+1] = AT_colptr
    sqrtδ = sqrt(δ)
    @inbounds for j=1:ncon
        @simd for k=AT_colptr[j]:(AT_colptr[j+1]-1)
            K_rowval[k] = AT_rowval[k]
            K_nzval[k] = AT_nzval[k] / sqrtδ
        end
    end

    # 2nd bloc
    K_deb = K_colptr[ncon+1]
    @inbounds @simd for i=1:nvar
        K_colptr[ncon+1+i] = K_deb + LH_colptr[i+1] - 1
    end
    @inbounds for j=1:nvar 
        @simd for k=LH_colptr[j]: (LH_colptr[j+1]-1)
            idx_K = K_deb + k - 1
            K_rowval[idx_K] = LH_rowval[k]
            K_nzval[idx_K] = LH_nzval[k]
        end
    end

    # third bloc
    K_deb = K_colptr[ncon+nvar+1]
    @inbounds @simd for i=1:nvar
        idx_K = K_deb + i 
        K_colptr[ncon+nvar+i+1] = idx_K
        K_rowval[idx_K-1] = i 
        K_nzval[idx_K-1] = sqrt(D[i])
    end
end

function fill_K1_2lp!(K_colptr, K_rowval, K_nzval, AT_colptr, AT_rowval, AT_nzval, D, δ, nvar, ncon)

    # fill first bloc
    K_colptr[1:ncon+1] = AT_colptr
    sqrtδ = sqrt(δ)
    @inbounds for j=1:ncon
        @simd for k=AT_colptr[j]:(AT_colptr[j+1]-1)
            K_rowval[k] = AT_rowval[k]
            K_nzval[k] = AT_nzval[k] / sqrtδ
        end
    end

    # 2nd bloc
    K_deb = K_colptr[ncon+1]
    @inbounds @simd for i=1:nvar
        idx_K = K_deb + i 
        K_colptr[ncon+i+1] = idx_K
        K_rowval[idx_K-1] = i 
        K_nzval[idx_K-1] = sqrt(D[i])
    end
end

function update_K1_2!(K_nzval, AT_colptr, AT_nzval, D, δ, ncon, nvar)

    sqrtδ = sqrt(δ)
    @inbounds for j=1:ncon
        @simd for k=AT_colptr[j]:(AT_colptr[j+1]-1)
            K_nzval[k] = AT_nzval[k] / sqrtδ
        end
    end

    K_nzval[end-nvar+1:end] .= sqrt.(D)
end