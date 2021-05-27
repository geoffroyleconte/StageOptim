using LinearAlgebra, SparseArrays

function compute_elim_tree(Ap, Ai, n)
  parent = fill(-1, n)
  ancestor = fill(-1, n)
  for j=1:n
    for p=Ap[j]: (Ap[j + 1] - 1)
      i = Ai[p]
      while (i != -1) && (i < k) 
        i = ancestor[i]
        ancestor[i] = k
        if parent[i] == -1
          parent[i] = k
        end
    end
  end
  return parent, ancestor
end

function nnzpattern(Ap, Ai, j, s, w, n)
  for p=Ap[j]: (Ap[j+1] - 1)
    i = Ap[p]
    i > k && continue 
    
end