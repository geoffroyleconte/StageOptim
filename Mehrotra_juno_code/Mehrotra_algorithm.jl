using LinearAlgebra
using DataStructures
using QuadraticModels
using SparseArrays
using NLPModels
using Printf

function display_results(result)
    # fonction pour l'affichage
    println("\n------------------------------------------------------------------------------------------------------")
    println("---------------------------------------------- RESULTS -----------------------------------------------")
    display(result)
end

function Compute_AlphaAff_dual(alpha_step, v_k, dir_v_k)
    alpha = 0
    n = length(v_k)
    while alpha+alpha_step <= 1 && all(v_k + (alpha+alpha_step) * dir_v_k .>= 0)
        alpha += alpha_step
    end
    return alpha
end

function Compute_AlphaMax_dual(alpha_step, v_k, dir_v_k)
    alpha = 0
    n = length(v_k)
    while alpha+alpha_step <= 100/99 && all(v_k + (alpha+alpha_step) * dir_v_k .>= 0)
        alpha += alpha_step
    end
    return alpha
end

function Compute_AlphaAff_primal(alpha_step, v_k, dir_v_k, lvar, uvar)
    alpha = 0
    n = length(v_k)
    while alpha+alpha_step <= 1 && all(lvar .<= (v_k + (alpha+alpha_step) * dir_v_k) .<= uvar)
        alpha += alpha_step
    end
    return alpha
end

function Compute_AlphaMax_primal(alpha_step, v_k, dir_v_k, lvar, uvar)
    alpha = 0
    n = length(v_k)
    while alpha+alpha_step <= 100/99 && all(lvar .<= (v_k + (alpha+alpha_step) * dir_v_k) .<= uvar)
        alpha += alpha_step
    end
    return alpha
end

function Compute_mu(x_l, x_u, s, stilde, lvar, uvar)
    #x_l coordinates of x corresponding to finite lower bounds ( resp. finite upper bounds for x_u)
    # arguments must have finite bounds
    return (s' * x_l + stilde' * x_u) / (length(x_l) + length(x_u))
end


function MehrotraPCQuadBounds(QM; max_iter=20, eps=1e-10, tol_step_x=1e-8, eps_mu=1e-8, alpha_step=1e-3, display=true)

    # get variables from QuadraticModel
    x_0 = QM.meta.x0 .- 0.5 # to have a starting point in [lvar, uvar]
    #x_0 = vcat(ones(32)*1, ones(19)*(-1))
    lvar, uvar = QM.meta.lvar, QM.meta.uvar
    id_non_inf_lvar, id_non_inf_uvar = findall((x -> x!=-Inf), lvar), findall((x -> x!=Inf), uvar) # finite bounds index
    @assert all(x_0 .> lvar) && all(x_0 .< uvar)
    A = jac(QM, x_0)
    n_rows, n_cols = size(A)
    Q = hess(QM, x_0)
    c = grad(QM, zeros(n_cols))
    c0 = obj(QM, zeros(n_cols))
    @assert QM.meta.lcon == QM.meta.ucon # equality constraint (Ax=b)
    b = QM.meta.lcon
    s_0, stilde_0 = zeros(n_cols), zeros(n_cols)
    s_0[id_non_inf_lvar] = ones(length(id_non_inf_lvar))
    stilde_0[id_non_inf_uvar] = ones(length(id_non_inf_uvar))
    lambda_0 = Matrix(A)'\(c -s_0+stilde_0) # least square initialisation, s_0 = stilde_0
    x_k, lambda_k, s_k, stilde_k = copy(x_0), copy(lambda_0), copy(s_0), copy(stilde_0)

    rb_0 = A * x_0 - b
    rc_0 = -Q * x_0 + A' * lambda_0 + s_0 - stilde_0 - c
    mu_0 = Compute_mu(x_0[id_non_inf_lvar], x_0[id_non_inf_uvar],
                        s_0[id_non_inf_lvar], stilde_0[id_non_inf_uvar],
                        lvar[id_non_inf_lvar], uvar[id_non_inf_uvar])

    # matrices without infinity constraints
    nb_non_inf_l, nb_non_inf_u = length(id_non_inf_lvar), length(id_non_inf_uvar) # number of finite constraints
    nb_inf_l, nb_inf_u = n_cols-nb_non_inf_l, n_cols-nb_non_inf_u # number of infinite constraints
    Xk_non_inf_l = Diagonal(x_k[id_non_inf_lvar])
    Sk_non_inf = [Diagonal(s_k[id_non_inf_lvar])  zeros(nb_non_inf_l, nb_inf_l)] # add zeros to match number of cols
    Xk_non_inf_u = Diagonal(x_k[id_non_inf_uvar])
    Stildek_non_inf = [Diagonal(stilde_k[id_non_inf_uvar])  zeros(nb_non_inf_u, nb_inf_u)]
    Lvar_non_inf = Diagonal(lvar[id_non_inf_lvar])
    Uvar_non_inf = Diagonal(uvar[id_non_inf_uvar])

    mu_k, rb_k, rc_k = copy(mu_0), copy(rb_0), copy(rc_0)
    k = 0
    e = ones(n_cols)
    n_c = norm(c)
    n_b = norm(b)

    # stopping criterion
    quad_part = x_k' * Q * x_k
    pdd = abs(quad_part + c' * x_k - b' * lambda_k ) / (1 + abs(c' * x_k + quad_part/2))
    cond_rb, cond_rc = norm(rb_k) / (1 + n_b), norm(rc_k) / (1 + n_c)
    opti_pdd, opti_rb, opti_rc = pdd < eps, cond_rb < eps, cond_rc < eps
    small_step_x, small_mu = false, mu_k < eps_mu

    # display
    if display == true
        println("Iter | primal_objective | primal-dual difference | rb condition | rc condition |   step x   |     mu")
        println("--------------------------------------------------------------------------------------------------------")
        @printf("% 4d |     % 7.2e    |        % 7.2e       |   % 7.2e  |   % 7.2e  | % 7.2e  | % 7.2e\n",
                k, c0+c'*x_k +quad_part/2, pdd, cond_rb, cond_rc,0., mu_k)
    end


    while k<max_iter && opti_pdd==false && opti_rb==false && opti_rc==false && small_step_x==false && small_mu==false

        # Affine scaling direction
        Jacob_Fk = [-Q                      A'                     I(n_cols)[:, 1:nb_non_inf_l]       I(n_cols)[:, 1:nb_non_inf_u]*(-1)
                    A                zeros(n_rows, n_rows)         zeros(n_rows, nb_non_inf_l)        zeros(n_rows, nb_non_inf_u)
                    Sk_non_inf       zeros(nb_non_inf_l, n_rows)    Xk_non_inf_l-Lvar_non_inf         zeros(nb_non_inf_l, nb_non_inf_u)
                    Stildek_non_inf  zeros(nb_non_inf_u, n_rows)   zeros(nb_non_inf_u,nb_non_inf_l)   Xk_non_inf_u-Uvar_non_inf]

        Fk_aff = [-rc_k
                  -rb_k
                  -(x_k[id_non_inf_lvar]-lvar[id_non_inf_lvar]).*s_k[id_non_inf_lvar]
                  -(x_k[id_non_inf_uvar]-uvar[id_non_inf_uvar]).*stilde_k[id_non_inf_uvar]]

        dir_aff_k = Jacob_Fk\Fk_aff

        alpha_aff_pri = Compute_AlphaAff_primal(alpha_step, x_k, dir_aff_k[1:n_cols], lvar, uvar)
        alpha_aff_dual = Compute_AlphaAff_dual(alpha_step, s_k[id_non_inf_lvar],
                                            dir_aff_k[n_rows+n_cols+1: n_rows+n_cols+nb_non_inf_l])
        alphatilde_aff_dual = Compute_AlphaAff_dual(alpha_step, stilde_k[id_non_inf_uvar],
                                                dir_aff_k[n_rows+n_cols+nb_non_inf_l+1:end])

        mu_aff = Compute_mu(x_k[id_non_inf_lvar] + alpha_aff_pri * dir_aff_k[1:n_cols][id_non_inf_lvar],
                    x_k[id_non_inf_uvar] + alpha_aff_pri * dir_aff_k[1:n_cols][id_non_inf_uvar],
                    s_k[id_non_inf_lvar] + alpha_aff_dual * dir_aff_k[n_rows+n_cols+1: n_rows+n_cols+nb_non_inf_l],
                    stilde_k[id_non_inf_uvar] + alphatilde_aff_dual * dir_aff_k[n_rows+n_cols+nb_non_inf_l+1: end],
                    lvar[id_non_inf_lvar], uvar[id_non_inf_uvar])

        sigma = (mu_aff / mu_k)^3

        # corrector and centering step
        Fk_cc = [zeros(n_rows+n_cols, 1)
                 sigma*mu_k*e[1:nb_non_inf_l] - dir_aff_k[1:n_cols][id_non_inf_lvar].*dir_aff_k[n_rows+n_cols+1: n_rows+n_cols+nb_non_inf_l]
                 sigma*mu_k*e[1:nb_non_inf_u] - dir_aff_k[1:n_cols][id_non_inf_uvar].*dir_aff_k[n_rows+n_cols+nb_non_inf_l+1: end]]
        dir_cc_k = Jacob_Fk\Fk_cc

        dir_k = dir_aff_k .+ dir_cc_k # final direction

        alpha_max_pri = Compute_AlphaMax_primal(alpha_step, x_k, dir_k[1:n_cols], lvar, uvar)
        alpha_max_dual = Compute_AlphaMax_dual(alpha_step, s_k[id_non_inf_lvar],
                                        dir_k[n_rows+n_cols+1: n_rows+n_cols+nb_non_inf_l])

        alphatilde_max_dual = Compute_AlphaMax_dual(alpha_step, stilde_k[id_non_inf_uvar],
                                            dir_k[n_rows+n_cols+nb_non_inf_l+1: end])

        # new parameters
        alpha_k_pri = min(0.99*alpha_max_pri, 1)
        alpha_k_dual = min(0.99*alpha_max_dual, 1)
        alphatilde_k_dual = min(0.99*alphatilde_max_dual, 1)
        x_k += alpha_k_pri * dir_k[1:n_cols]
        lambda_k += alpha_k_dual * dir_k[n_cols+1: n_rows+n_cols]
        s_k[id_non_inf_lvar] += alpha_k_dual * dir_k[n_rows+n_cols+1: n_rows+n_cols+nb_non_inf_l]
        stilde_k[id_non_inf_uvar] += alphatilde_k_dual * dir_k[n_rows+n_cols+nb_non_inf_l+1: end]

        Xk_non_inf_l = Diagonal(x_k[id_non_inf_lvar])
        Sk_non_inf = [Diagonal(s_k[id_non_inf_lvar])  zeros(nb_non_inf_l, nb_inf_l)]
        Xk_non_inf_u = Diagonal(x_k[id_non_inf_uvar])
        Stildek_non_inf = [Diagonal(stilde_k[id_non_inf_uvar])  zeros(nb_non_inf_u, nb_inf_u)]

        step_x =  10 #norm(alpha_k_pri * dir_k[1:n_cols])
        mu_k = Compute_mu(x_k[id_non_inf_lvar], x_k[id_non_inf_uvar],
                        s_k[id_non_inf_lvar], stilde_k[id_non_inf_uvar],
                        lvar[id_non_inf_lvar], uvar[id_non_inf_uvar])

        rb_k = A * x_k - b
        rc_k = -Q * x_k + A' * lambda_k + s_k - stilde_k - c

        # update stopping criterion values:
        quad_part = x_k' * Q * x_k
        pdd = abs(quad_part + c' * x_k - b' * lambda_k ) / (1 + abs(c' * x_k + quad_part/2)) # test correct? ecart primal dual
        cond_rb = norm(rb_k) / (1 + n_b)
        cond_rc = norm(rc_k) / (1 + n_c)
        opti_pdd, opti_rb, opti_rc = pdd < eps, cond_rb < eps, cond_rc < eps
        small_step_x, small_mu = step_x < tol_step_x, mu_k < eps_mu

        k += 1

        if display == true
            @printf("% 4d |     % 7.2e    |        % 7.2e       |   % 7.2e  |   % 7.2e  | % 7.2e  | % 7.2e\n",
                k, c0+c'*x_k +quad_part/2, pdd, cond_rb, cond_rc,step_x, mu_k)
        end
    end

    if display == true
        criteria = [k >= max_iter,  opti_pdd, opti_rb, opti_rc, small_step_x, small_mu]
        criteria_names = ["reached max_iter",  "pdd <= eps", "cond_rb <= eps", "cond_rc <= eps",
            "step_x <= small_step_x", "mu_k <= eps_mu"]
        println("\n stopping criterion = ",criteria_names[findall(criteria)])
    end

    return OrderedDict("x_opt" => x_k, "lambda_opt" => lambda_k, "s_opt" => s_k, "stilde_opt" => stilde_k,
        "n_iter" => k, "pdd" => pdd, "cond_rb" => cond_rb, "cond_rc" => cond_rc)
end
