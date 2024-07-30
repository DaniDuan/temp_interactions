"
To include temperature function 
"

function temp_trait(N, kw) 
    k = 0.0000862 # Boltzman constant
    @unpack T, Tr, Ed= kw
    B,E,Tp = randtemp_param(N, kw)
    temp_p = B .* exp.((-E./k) * ((1/T)-(1/Tr)))./(1 .+ (E./(Ed .- E)) .* exp.(Ed/k * (1 ./Tp .- 1/T)))
    # Eϵ = (B[:,2] .* (E[:,1] .- E[:,2]))./(B[:,1]*(1-L)-B[:,2])  #m0(Eu − Em)/(u0(1 − l) − m0)
    return temp_p, E, Tp
end

function randtemp_param(N, kw)
    @unpack T, ρ_t, Tr, Ed = kw
    k = 0.0000862 # Boltzman constant
    # Here setting B0_u = m0 /(1 - L - CUE_0) = 0.2772
    # with L = 0.3, mean(m0) = 0.132, CUE_0 = 0.22
    B0 = [-0.8116 -1.4954]# Using CUE0 = 0.22, mean growth rate = 0.48
    B0_var = 0.17 .* abs.(B0); E_mean = [0.8146 0.5741]; E_var =  0.1364 .* E_mean
    cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
    meanv = [B0 ; E_mean]
    cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
    cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]
    allu = rand(MvNormal(meanv[:,1], cov_u), N)
    allm = rand(MvNormal(meanv[:,2], cov_m), N)
    B = [exp.(allu[1,:]) exp.(allm[1,:])]
    E = [allu[2,:] allm[2,:]]
    
    Tpu = 273.15 .+ rand(Normal(35, 5), N)
    Tpm = Tpu .+ 3
    Tp = [Tpu Tpm]
    return B,E,Tp
end 
