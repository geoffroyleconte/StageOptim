using LinearOperators, LimitedLDLFactorization

function lldl_op(K :: SparseMatrixCSC, p :: Int; droptol :: Float64=0.0, minpiv :: Float64=0.0)
  n = size(K, 1);
  (Lmat, d, alpha) = lldl(K, p, droptol=droptol, minpiv=minpiv);
  D = opDiagonal(1./abs(d));
  L = opInverse(Lmat + speye(n));
  return (L' * D * L, alpha);
end