
function Eff_Lv_Jac(; p_lv, sol)

    N = p_lv.N
    ℵ = p_lv.ℵ
    r = p_lv.r
    C = sol[1:N, length(sol)]
    # Calculating Jacobian matrix
    LV_Jac = [ℵ[i, j]*C[i] for i in 1:N, j in 1:N]
    # reset diagonal
    LV_Jac[diagind(LV_Jac)] .= [r[i] + 2*ℵ[i, i]*C[i] + sum(ℵ[i, j]*C[j] for j in 1:N if j != i) for i in 1:N]

    return(LV_Jac)
end
