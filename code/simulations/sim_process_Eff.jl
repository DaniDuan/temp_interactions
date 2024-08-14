include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob
path = glob("Eff_iters*", "../data/Eff_p0/")

N=100
M=50
### Temp params 
# ρ_t= [-0.3500 -0.3500]; # realistic covariance
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

Temp_rich = range(0, num_temps-1, length = num_temps)

progress = Progress(length(path)*num_temps; desc="Progress running:")

all_Rrela_collect = Vector{Vector{Float64}}(); all_Crela_collect = Vector{Vector{Float64}}()
Eff_results = zeros(Float64, num_temps, 59)
@time for j in 1: num_temps
    all_ℵii_H = Float64[]; all_ℵij_H = Union{Float64, Missing}[]; all_ℵij_d_H = Union{Float64, Missing}[];
    all_uℵij_H = Union{Float64, Missing}[]; all_lℵij_H = Union{Float64, Missing}[]; 
    all_ℵii_sur_H = Union{Float64, Missing}[]; all_ℵij_sur_H = Union{Float64, Missing}[]; all_ℵij_d_sur_H = Union{Float64, Missing}[];
    all_uℵij_sur_H = Union{Float64, Missing}[]; all_lℵij_sur_H = Union{Float64, Missing}[];
    all_r_H = Float64[]; 
    all_leading_H = Float64[]; all_diag_H = Float64[];radi_H = Float64[]; diag_dominance_H = Float64[];
    all_u_H = Float64[]; all_m_H = Float64[]; RO_H = Union{Float64, Missing}[]; ulO_H = Union{Float64, Missing}[]; Rul_H = Union{Float64, Missing}[]; 
    all_Eu_H = Float64[]; all_Em_H = Float64[]; all_Eu_sur_H = Float64[]; all_Em_sur_H = Float64[];
    all_Tpu_H = Float64[]; all_Tpm_H = Float64[]; all_Tpu_sur_H = Float64[]; all_Tpm_sur_H = Float64[]; all_sumαij = Float64[];
    all_Rrela_H = Float64[]; all_Crela_H = Float64[]

    for i in 1:length(path)
        @load path[i] all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_ℵii_sur all_ℵij_sur all_ℵij_d_sur all_uℵij_sur all_lℵij_sur all_r  all_u all_m RO ulO Rul all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_diag radi diag_dominance all_Rrela all_Crela
        append!(all_ℵii_H, all_ℵii[j]); append!(all_ℵij_H, all_ℵij[j]); append!(all_ℵij_d_H,all_ℵij_d[j]); append!(all_uℵij_H,all_uℵij[j]); append!(all_lℵij_H,all_lℵij[j]);
        append!(all_ℵii_sur_H, all_ℵii_sur[j]); append!(all_ℵij_sur_H, all_ℵij_sur[j]); append!(all_ℵij_d_sur_H, all_ℵij_d_sur[j]); append!(all_uℵij_sur_H, all_uℵij_sur[j]); append!(all_lℵij_sur_H, all_lℵij_sur[j]);
        append!(all_r_H, all_r[j]);
        append!(all_u_H, all_u[j]); append!(all_m_H, all_m[j]); append!(RO_H, RO[j]); append!(ulO_H, ulO[j]); append!(Rul_H, Rul[j]);
        append!(all_Eu_H, all_Eu[j]); append!(all_Em_H, all_Em[j]); append!(all_Eu_sur_H, all_Eu_sur[j]); append!(all_Em_sur_H, all_Em_sur[j]);
        append!(all_Tpu_H, all_Tpu[j]); append!(all_Tpm_H, all_Tpm[j]); append!(all_Tpu_sur_H, all_Tpu_sur[j]); append!(all_Tpm_sur_H, all_Tpm_sur[j]);
        push!(all_leading_H, all_leading[j]); append!(all_diag_H, all_diag[j]); append!(radi_H, radi[j]); push!(diag_dominance_H, diag_dominance[j]);
        append!(all_Rrela_H, all_Rrela[j]); append!(all_Crela_H, all_Crela[j])    

        next!(progress)

    end 
    push!(all_Rrela_collect, all_Rrela_H); push!(all_Crela_collect,all_Crela_H)
    Eff_results[Int(j),:] = [mean(all_ℵii_H), std(all_ℵii_H)/sqrt(length(all_ℵii_H)), mean(skipmissing(all_ℵij_H)), std(skipmissing(all_ℵij_H))/sqrt(length(all_ℵij_H)), 
    mean(skipmissing(all_ℵij_d_H)), std(skipmissing(all_ℵij_d_H))/sqrt(length(all_ℵij_d_H)), mean(skipmissing(all_uℵij_H)), std(skipmissing(all_uℵij_H))/sqrt(length(all_uℵij_H)), mean(skipmissing(all_lℵij_H)), std(skipmissing(all_lℵij_H))/sqrt(length(all_lℵij_H)), 
    mean(all_ℵii_sur_H), std(all_ℵii_sur_H)/sqrt(length(all_ℵii_sur_H)), mean(skipmissing(all_ℵij_sur_H)), std(skipmissing(all_ℵij_sur_H))/sqrt(length(all_ℵij_sur_H)), 
    mean(skipmissing(all_ℵij_d_sur_H)), std(skipmissing(all_ℵij_d_sur_H))/sqrt(length(all_ℵij_d_sur_H)), mean(skipmissing(all_uℵij_sur_H)), std(skipmissing(all_uℵij_sur_H))/sqrt(length(all_uℵij_sur_H)), mean(skipmissing(all_lℵij_sur_H)), std(skipmissing(all_lℵij_sur_H))/sqrt(length(all_lℵij_sur_H)), 
    mean(all_r_H), std(all_r_H)/sqrt(length(all_r_H)), mean(all_u_H), std(all_u_H)/sqrt(length(all_u_H)), mean(all_m_H), std(all_m_H)/sqrt(length(all_m_H)), 
    mean(skipmissing(RO_H)), std(skipmissing(RO_H))/sqrt(length(RO_H)),
    mean(skipmissing(ulO_H)), std(skipmissing(ulO_H))/sqrt(length(ulO_H)),
    mean(skipmissing(Rul_H)), std(skipmissing(Rul_H))/sqrt(length(Rul_H)),
    mean(all_Eu_H), std(all_Eu_H)/sqrt(length(all_Eu_H)), mean(all_Em_H), std(all_Em_H)/sqrt(length(all_Em_H)),
    mean(all_Eu_sur_H), std(all_Eu_sur_H)/sqrt(length(all_Eu_sur_H)), mean(all_Em_sur_H), std(all_Em_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_Tpu_H), std(all_Tpu_H)/sqrt(length(all_Tpu_H)), mean(all_Tpm_H), std(all_Tpm_H)/sqrt(length(all_Tpm_H)),
    mean(all_Tpu_sur_H), std(all_Tpu_sur_H)/sqrt(length(all_Tpm_sur_H)), mean(all_Tpm_sur_H), std(all_Tpm_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_leading_H), std(all_leading_H)/sqrt(length(all_leading_H)), sum(real.(all_leading_H) .< 0)/length(path), 
    mean(diag_dominance_H), std(diag_dominance_H)/sqrt(length(diag_dominance_H)), 
    mean(all_diag_H), std(all_diag_H)/sqrt(length(all_diag_H)), mean(radi_H), std(radi_H)/sqrt(length(radi_H)),
    mean(all_sumαij), std(all_sumαij)/sqrt(length(all_sumαij))]
end 
R"library(beepr); beep(sound = 4, expr = NULL)"
col_names_EF = ["αii", "αii_err", "αij", "αij_err", "αij_d", "αij_d_err", "αij_upper", "αij_upper_err","αij_lower", "αij_lower_err",
                "αii_sur", "αii_sur_err", "αij_sur", "αij_sur_err", "αij_d_sur", "αij_d_sur_err", "αij_upper_sur", "αij_upper_sur_err","αij_lower_sur", "αij_lower_sur_err",
                "r", "r_err", "u", "u_err","m", "m_err", 
                "RO", "RO_err", "ulO", "ulO_err", "estα", "estα_err",
                "Eu", "Eu_err", "Em", "Em_err", "Eu_sur", "Eu_sur_err", "Em_sur", "Em_sur_err",
                "Tpu", "Tpu_err", "Tpm", "Tpm_err", "Tpu_sur", "Tpu_sur_err", "Tpm_sur", "Tpm_sur_err",
                "eigen", "eigen_err", "stability",
                "diag_dom", "diag_dom_err",
                "Jac_diag", "Jac_diag_err", "radius", "radius_err", "sum_αij", "sum_αij_err"];
Eff_results = DataFrame(Eff_results, col_names_EF);

# CSV.write("../results/Eff_results_p0.csv", Eff_results, writeheader=false)
# @save "../results/Feas_CR_dist_p0.jld2" all_Rrela_collect all_Crela_collect
# @load "../results/Feas_CR_dist_p0.jld2" all_Rrela_collect all_Crela_collect


temp = collect(Temp_rich .+273.15)

all_temp_R = vcat([repeat([temp[t]], length(all_Rrela_collect[t])) for t in 1:num_temps]...)
[mean(all_Rrela_collect[t]) for t in 1:num_temps]
f = Figure(resolution = (1200, 1200));
ax1 = Axis(f[1,1], xlabel = "Temperature", ylabel = "Resource distribution", ygridvisible = false, xgridvisible = false)
scatter!(ax1,all_temp_R .- 273.15, vcat(all_Rrela_collect...), color = ("#285C93"), label = "", alpha = 0.7)

f

all_temp_C = vcat([repeat([temp[t]], length(all_Crela_collect[t])) for t in 1:num_temps]...)
f = Figure(resolution = (1200, 1200));
ax2 = Axis(f[1,1], xlabel = "Temperature", ylabel = "Consumer distribution", ygridvisible = false, xgridvisible = false)
scatter!(ax2,all_temp_C .- 273.15, vcat(all_Crela_collect...), color = ("#285C93"), label = "", alpha = 0.7)
f


save("../results/Feas_R_dist_p0.png", f) 

f = Figure(resolution = (1200, 1200));
for j in 1: num_temps
    # if (j-1) % 5 == 0
        circ_leading = Float64[]
        for i in 1:length(path)
            @load path[i] all_leading
            push!(circ_leading, all_leading[j])
        end 
        ax = Axis(f[Int(floor((j-1)/5+1)),Int((j-1) % 5+1)], xlabel = "real", ylabel = "imaginary", title = "T = $(j-1) °C", ygridvisible = false, xgridvisible = false)
        xlims!(ax, -0.25, 0.25)
        scatter!(ax, real.(circ_leading), imag.(circ_leading), color = ("#285C93"), label = "", markersize = 7, alpha = 0.7)
        lines!(ax, [0,0], [-1,1], linestyle = :dash, color = ("#4F363E", 1))
        # axislegend(position = :rb)
    # end 
end 
f
save("../results/leading_ρ-1.png", f) 


########## TPC of α ###########
all_Toii = Float64[]; all_Toij = Float64[]
for i in 1:length(path)
    @load path[i] all_ℵ all_ℵii
    vαii = [var(log.(abs.(all_ℵii[i]))) for i in 1:num_temps]
    Toii = Temp_rich[argmin(vαii)]
    αij = [sum(reshape(all_ℵ[i],Int(sqrt(length(all_ℵ[i]))),Int(sqrt(length(all_ℵ[i])))), dims = 2) for i in 1:31]
    vαij = [var(log.(abs.(αij[i]))) for i in 1:num_temps]
    Toij = Temp_rich[argmin(vαij)]
    push!(all_Toii, Toii); push!(all_Toij, Toij)
end 

Tt = Int(floor((mean(all_Toii) + mean(all_Toij))/2))

k = 0.0000862 # Boltzman constant
x = -1/k .* (1 ./(range(0, Tt, length = Tt+1) .+273.15) .- 1/Tr)
x_t = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)

f1 = Figure(fontsize = 35, resolution = (1200, 900));
f2 = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f1[1,1], xlabel = "Temperature (°C)", ylabel = "αii", xlabelsize = 50, ylabelsize = 50)
ax2 = Axis(f2[1,1], xlabel = "Temperature (°C)", ylabel = "αij", xlabelsize = 50, ylabelsize = 50)

all_Bαii = Float64[]; all_Bαij = Float64[]
all_Eαii = Float64[]; all_Eαij = Float64[]
for i in 1:length(path)
    @load path[i] all_ℵ all_ℵii

    y_αii = [mean(all_ℵii[i]) for i in 1:(Tt+1)]
    data_ii = DataFrame(y = log.(abs.(y_αii)), x = x);
    Bαii, Eαii= coef(lm(@formula(y ~ x), data_ii))
    Bαii = exp(Bαii)

    αij = [sum(reshape(all_ℵ[i],Int(sqrt(length(all_ℵ[i]))),Int(sqrt(length(all_ℵ[i])))), dims = 2).- diag(reshape(all_ℵ[i],Int(sqrt(length(all_ℵ[i]))),Int(sqrt(length(all_ℵ[i]))))) for i in 1:num_temps]
    all_yij = [mean(αij[i]) for i in 1:num_temps]
    y_αij = [mean(αij[i]) for i in 1:(Tt+1)]
    data_ij = DataFrame(y = log.(abs.(y_αij)), x = x);
    Bαij, Eαij = coef(lm(@formula(y ~ x), data_ij))
    Bαij = exp(Bαij)
    yii = Bαii * exp.(Eαii .* x_t)
    yij = Bαij * exp.(Eαij .* x_t)
    lines!(ax1, Temp_rich, log.(yii), color = ("#E17542", 0.2), linewidth = 0.5)
    lines!(ax2, Temp_rich, log.(yij), color = ("#E17542", 0.2), linewidth = 0.5)
    push!(all_Bαii, Bαii); push!(all_Bαij, Bαij)
    push!(all_Eαii, Eαii); push!(all_Eαij, Eαij)
end 

mean_ii = DataFrame(y = log.(abs.(Eff_results.αii)), x = x_t);
Bii, Eii= coef(lm(@formula(y ~ x), mean_ii))
Bii = exp(Bii)
yii = log.(Bii * exp.(E .* x_t))
lines!(ax1, Temp_rich, log.(abs.(Eff_results.αii)), color = ("#285C93", 1), linewidth = 7, label = "αii")
band!(ax1, Temp_rich, log.(abs.(Eff_results.αii .- Eff_results.αii_err)), log.(abs.(Eff_results.αii .+ Eff_results.αii_err)), color = ("#285C93", 0.5))
lines!(ax1, Temp_rich, yii, color = ("#4F363E", 1), linewidth = 5, label = "")

mean_ij = DataFrame(y = log.(abs.(Eff_results.sum_αij))[1:Tt+1], x = x);
Bij, Eij= coef(lm(@formula(y ~ x), mean_ij))
Bij = exp(Bij)
yij = log.(Bij * exp.(Eij .* x_t))
lines!(ax2, Temp_rich, log.(abs.(Eff_results.sum_αij)), color = ("#285C93", 1), linewidth = 7, label = "αij")
band!(ax2, Temp_rich,  log.(abs.(Eff_results.sum_αij .- Eff_results.sum_αij_err)), log.(abs.(Eff_results.sum_αij .+ Eff_results.sum_αij_err)), color = ("#285C93", 0.5))
lines!(ax2, Temp_rich, yij, color = ("#4F363E", 1), linewidth = 5, label = "")

f1
save("../results/temp_αii.png", f1) 

f2
save("../results/temp_αij.png", f2) 
