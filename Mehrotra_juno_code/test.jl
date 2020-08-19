include("Mehrotra_algorithm.jl")

# tests pour le presolve
A = sparse([0. 0. 0.;
            0. 0. 0.;
            2. 2. 2. ;
            0. 0. 0.;
            0. 0. 0.;
            0. 1. 0.;
            0. 2. 2.;
            0. 0. 0.;
            1. 1. 1.;
            0. 0. 0.;
            0. 0. 0.])
b = [0.; 0. ;1. ; 0.; 0.; 0.; 2.; 0.; 2.; 0.; 0.]

Q = sparse([1. 0. 0.;
            2. 3. 0.;
            0. 2. 4.])
c = [1. ;2.; 3.]
c0 = 1.
lvar = [-Inf; 2.; -Inf]
uvar = [Inf; 20; 18.]
ifix = []
ilow = findall(x->x!=-Inf, lvar)
iupp = findall(x->x!=Inf, uvar)
irng = [2]
dropzeros!(A)
dropzeros!(Q)
Arows, Acols, Avals = findnz(A)
Qrows, Qcols, Qvals = findnz(Q)

n_rows = size(A)[1]
n_cols = size(A)[2]
Arows, Acols, Avals, b, lvar, uvar, ifix = rm_singleton_rows!(Arows, Acols, Avals,
                                                              b, lvar, uvar, ifix, n_rows)
e_r = find_empty_rows(Arows, n_rows)
Arows, b, n_rows = rm_empty_rows!(e_r, Arows, b, n_rows)

@code_warntype rm_ifix!(ifix, ilow, iupp, irng, Qrows, Qcols, Qvals, c, c0,
                                                  Arows, Acols, Avals, b, lvar, uvar, n_rows, n_cols)
Qrows, Qcols, Qvals, c, c0,
    Arows, Acols, Avals, b,
    lvar, uvar,ilow, iupp, irng, n_cols = rm_ifix!(ifix, ilow, iupp, irng, Qrows, Qcols, Qvals, c, c0,
                                                  Arows, Acols, Avals, b,lvar, uvar, n_rows, n_cols)


display(Matrix(sparse(Arows, Acols, Avals)))
display(Matrix(Symmetric(sparse(Qrows, Qcols, Qvals),:L)))
println(b, n_rows, n_cols)
println(lvar, uvar, ifix, c0, b)
println(ilow, iupp, irng)
