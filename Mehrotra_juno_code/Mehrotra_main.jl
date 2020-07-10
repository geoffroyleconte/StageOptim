include("Mehrotra_algorithm.jl")
#using ProfileView

# test with a quadratic problem
Q = [6. 2 1
    2 5 2
    1 2 4]
c = [-8.; -3; -3]
c0 = 0.
A = [1. 0 1
    0 1 1]
b = [0.; 3];
lvar = [0.;0;0]
uvar = [Inf; Inf; Inf]
lcon = b
ucon = b
x01 = [1.; 2.; 3.];

QM = QuadraticModel(c, Q, A=A, lcon=lcon, ucon=ucon, lvar=lvar, uvar=uvar, x0=x01, c0=c0, name="QM1")
stats_mpc1 =  mehrotraPCQuadBounds(QM)
println(stats_mpc1)

# path for netlib problems
path_pb = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes_netlib"

afiro = string(path_pb, "\\AFIRO.SIF")
qpdata2 = readqps(afiro)
QM2 = createQuadraticModel(qpdata2)
stats_mpc2 =  mehrotraPCQuadBounds(QM2)
println(stats_mpc2)


# problem 3   kb2    obj  -1.7499001299E+03
kb2 = string(path_pb, "\\KB2.SIF")
qpdata3 = readqps(kb2)
QM3 = createQuadraticModel(qpdata3)
res_mpc3 =  mehrotraPCQuadBounds(QM3);

# problem 4 SC50A  obj  -6.4575077059E+01
SC50A = string(path_pb, "\\CRE-B.SIF")
qpdata4 = readqps(SC50A)
QM4 = createQuadraticModel(qpdata4)
res_mpc4 =  mehrotraPCQuadBounds(QM4, max_iter=100);
