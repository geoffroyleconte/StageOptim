"""Scale the columns of a sparse matrix A by their Euclidean norm.
"""
function scale_ls!(A)
  s = ones(size(A, 2))
  for j = 1 : size(A, 2)
    i = A.colptr[j]
    k = A.colptr[j+1] - 1
    nj = i <= k ? norm(A.nzval[i:k]) : 0.0
    if nj > 0.0
      A.nzval[i:k] ./= nj
      s[j] = nj
    end
  end
  return s
end
