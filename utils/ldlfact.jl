using LinearAlgebra, SparseArrays

flip(i) = -i - 2
unflip(i) = (i < 0) ? flip(i) : i
marked(w, j) = w[j] < 0
function mark!(w, j)
  w[j] = flip(w[j])
end  

# sparse triangular G x = b , b = B[:,k]
function spreach(Gp, Gi, Bp, Bi, k, xi, pinv, n)
  top = n
  for p = Bp[k] : (Bp[k + 1] - 1)
    if !marked(Gp, Bi[p])
      top = dfs(Bi[p], G, top, xi, xi[n+1:end])
    end
  end
end

function dfs(j, Gp, Gi, top, xi, pstack, pinv)
  head  = 0
  while head >= 0
    j = xi[head]
    jnew = pinv[j] # ou j sans permutation
    if !marked(Gp, j)
      mark!(Gp, j) # Gp marqué
      pstack[head] = (jnew < 0) ? 0 : unflip(Gp[jnew])
    end
    done = 1 # node j done
    p2 = (jnew < 0) ? 0 : unflip(Gp[jnew])
    for p = pstack[head] : p2 - 1
      i = Gi[p]
      marked(Gp) && continue
      head += 1
      xi[head] = i
      done = 0
      break
    end
    if done 
      head -= 1
      top -= 1
      xi[top] = j
    end
  end
  return top
end

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


mutable struct CgSolver{T,S,S2} <: KrylovSolver{T,S}
  x  :: S
  r  :: S
  p  :: S
  Ap :: S
  z  :: S2
end

CgSolver{T,S}(x::S, r::S, p::S, Ap::S, z::S2) where {T,S,S2} = CgSolver{T,S,S2}(x, r, p, Ap, z)

function CgSolver{T,S}(n, m) where {T,S}
  x  = S(undef, n)
  r  = S(undef, n)
  p  = S(undef, n)
  Ap = S(undef, n)
  z  = nothing
  solver = CgSolver{T,S}(x, r, p, Ap, z)
  return solver
end