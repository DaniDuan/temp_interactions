
function Eff_Lv_Jac(; p_lv, sol)
    ### when ONLY considering feaible resultant community
    # bm = sol.u[length(sol.t)][1:p_lv.N]
    # sur = (1:p_lv.N)[bm .> 1.0e-7]
    # N = length(sur)
    # ℵ = [p_lv.ℵ[sur[i],sur[j]] for i in 1:N, j in 1:N]
    # r = p_lv.r[sur]
    # C = bm[sur]

    ### Collecting params
    N = p_lv.N
    ℵ = p_lv.ℵ
    r = p_lv.r
    C = sol[1:N, length(sol)]
    # Calculating Jacobian matrix
    LV_Jac = [ℵ[i, j]*C[i] for i in 1:N, j in 1:N]
    # reset diagonal
    LV_Jac[diagind(LV_Jac)] .= [r[i] + ℵ[i, i]*C[i] + sum(ℵ[i, j]*C[j] for j in 1:N) for i in 1:N]

    return(LV_Jac)
end
