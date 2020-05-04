include("Mehrotra_algorithm.jl")
include("QP_functions.jl")

# probleme1
Q = [6 2 1
    2 5 2
    1 2 4]
c = [-8; -3; -3]
c0 = 0.
A = [1 0 1
    0 1 1]
b = [0; 3];
lvar = [0;0;0]
uvar = [Inf; Inf; Inf]
lcon = b
ucon = b

x01 = [1.; 2.; 3.];

QM = QuadraticModel(c, Q, A=A, lcon=lcon, ucon=ucon, lvar=lvar, uvar=uvar, x0=x01, c0=c0, name="QM1")
SM = SlackModel(QM)

res_mpc1 =  MehrotraPCQuadBounds(SM, max_iter=20);
display_results(res_mpc1)

# probleme 2 afiro
path = "C:\\Users\\Geoffroy Leconte\\Documents\\cours\\TFE\\code\\problemes"
afiro = string(path, "\\AFIRO.SIF")

qpdata2 = readqps(afiro)
QM2 = createQuadraticModel(qpdata2, "QM2")
SM2 = SlackModel(QM2);
displayQuadraticModel(SM2)
x0_test = vcat(ones(32)*100, ones(19)*(-100))

res_mpc2 =  MehrotraPCQuadBounds(SM2, max_iter=20, alpha_step=1e-4);
display_results(res_mpc2)
