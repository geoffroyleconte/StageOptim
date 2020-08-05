include("Mehrotra_algorithm.jl")
#using ProfileView

# # test with a quadratic problem
# Q = [6. 2 1
#     2 5 2
#     1 2 4]
# c = [-8.; -3; -3]
# c0 = 0.
# A = [1. 0 1
#     0 1 1]
# b = [0.; 3];
# lvar = [0.;0;0]
# uvar = [Inf; Inf; Inf]
# lcon = b
# ucon = b
# x01 = [1.; 2.; 3.];

# QM = QuadraticModel(c, Q, A=A, lcon=lcon, ucon=ucon, lvar=lvar, uvar=uvar, x0=x01, c0=c0, name="QM1")
# stats_mpc1 =  mehrotraPCQuadBounds(QM)
# println(stats_mpc1)

# path for netlib problems
netlib_path = fetch_netlib()
mm_path = fetch_mm()

# afiro = string(joinpath(path_pb, "AFIRO.SIF"))
prob = string(joinpath(mm_path, "CVXQP3_L.SIF"))
qpdata = readqps(prob)
qp = createQuadraticModel(qpdata)
stats_mpc =  mehrotraPCQuadBounds(qp, max_time=Inf)
println(stats_mpc)


# # problem 3   kb2    obj  -1.7499001299E+03
# kb2 = string(path_pb, "\\KB2.SIF")
# qpdata3 = readqps(kb2)
# QM3 = createQuadraticModel(qpdata3)
# res_mpc3 =  mehrotraPCQuadBounds(QM3);

# # problem 4 SC50A  obj  -6.4575077059E+01
# SC50A = string(path_pb, "\\CRE-B.SIF")
# qpdata4 = readqps(SC50A)
# QM4 = createQuadraticModel(qpdata4)
# res_mpc4 =  mehrotraPCQuadBounds(QM4, max_iter=100);
