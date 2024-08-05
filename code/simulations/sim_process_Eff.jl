include("./sim_frame.jl")

using Glob
path = glob("Eff_iters*", "../data/output/")

N=100
M=50
### Temp params 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
Ci = fill(0.1, N)
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 31


Eff_results = zeros(Float64, num_temps, 45)
@time for j in 1: num_temps

    all_ℵ_H = Float64[]; all_ℵii_H = Float64[]; allsur_ℵij_H = Union{Float64, Missing}[]; all_r_H = Float64[]; 
    all_leading_H = Float64[]; all_diag_H = Float64[];radi_H = Float64[]; diag_dominance_H = Float64[];
    all_u_H = Float64[]; all_m_H = Float64[]; RO_H = Union{Float64, Missing}[]; ulO_H = Union{Float64, Missing}[]; Rul_H = Union{Float64, Missing}[]; all_UDLD_H = Union{Float64, Missing}[];
    all_Eu_H = Float64[]; all_Em_H = Float64[]; all_Eu_sur_H = Float64[]; all_Em_sur_H = Float64[];
    all_Tpu_H = Float64[]; all_Tpm_H = Float64[]; all_Tpu_sur_H = Float64[]; all_Tpm_sur_H = Float64[]

    for i in 1:length(path)
        @load path[i] all_ℵ all_ℵii all_ℵij all_r  all_u all_m RO ulO Rul all_UDLD all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_diag radi diag_dominance 
        append!(all_ℵ_H, sur_ℵ); append!(all_ℵii_H, ℵii); append!(allsur_ℵij_H, ℵij); append!(all_r_H, sur_r);
        append!(all_u_H, u); append!(all_m_H, m); append!(RO_H, R_over); append!(ulO_H, ul_over); append!(Rul_H, Rul_over); append!(all_UDLD_H, UDLD);
        append!(all_Eu_H, Eu); append!(all_Em_H, Em); append!(all_Eu_sur_H, Eu_sur); append!(all_Em_sur_H, Em_sur);
        append!(all_Tpu_H, Tpu); append!(all_Tpm_H, Tpm); append!(all_Tpu_sur_H, Tpu_sur); append!(all_Tpm_sur_H, Tpm_sur);
        push!(all_leading_H, leading); append!(all_diag_H, jac_diag); append!(radi_H, jac_off); push!(diag_dominance_H, diag_dom)
    end 
    Eff_results[Int(i+1),:] = [mean(all_ℵ_H), std(all_ℵ_H)/sqrt(length(all_ℵ_H)), 
    mean(all_ℵii_H), std(all_ℵii_H)/sqrt(length(all_ℵii_H)), mean(skipmissing(allsur_ℵij_H)), std(skipmissing(allsur_ℵij_H))/sqrt(length(allsur_ℵij_H)), 
    mean(all_r_H), std(all_r_H)/sqrt(length(all_r_H)), mean(all_u_H), std(all_u_H)/sqrt(length(all_u_H)), 
    mean(all_m_H), std(all_m_H)/sqrt(length(all_m_H)), 
    mean(skipmissing(RO_H)), std(skipmissing(RO_H))/sqrt(length(RO_H)),
    mean(skipmissing(ulO_H)), std(skipmissing(ulO_H))/sqrt(length(ulO_H)),
    mean(skipmissing(Rul_H)), std(skipmissing(Rul_H))/sqrt(length(Rul_H)),
    mean(skipmissing(all_UDLD_H)), std(skipmissing(all_UDLD_H))/sqrt(length(all_UDLD_H)),
    mean(all_Eu_H), std(all_Eu_H)/sqrt(length(all_Eu_H)), mean(all_Em_H), std(all_Em_H)/sqrt(length(all_Em_H)),
    mean(all_Eu_sur_H), std(all_Eu_sur_H)/sqrt(length(all_Eu_sur_H)), mean(all_Em_sur_H), std(all_Em_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_Tpu_H), std(all_Tpu_H)/sqrt(length(all_Tpu_H)), mean(all_Tpm_H), std(all_Tpm_H)/sqrt(length(all_Tpm_H)),
    mean(all_Tpu_sur_H), std(all_Tpu_sur_H)/sqrt(length(all_Tpm_sur_H)), mean(all_Tpm_sur_H), std(all_Tpm_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_leading_H), std(all_leading_H)/sqrt(length(all_leading_H)), sum(all_leading_H .< 0)/length(path), 
    mean(diag_dominance_H), std(diag_dominance_H)/sqrt(length(diag_dominance_H)), 
    mean(all_diag_H), std(all_diag_H)/sqrt(length(all_diag_H)), mean(radi_H), std(radi_H)/sqrt(length(radi_H))]
    
    print(j-1, " °C Complete, ", "α ",mean(all_ℵ_H),"\n") 
end 