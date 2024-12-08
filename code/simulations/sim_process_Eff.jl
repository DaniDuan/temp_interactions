include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob
using ColorSchemes
path = glob("Eff_iters*", "../data/Eff_p-1_new/")
# path = glob("Eff_iters*", "../data/L07/p-1/")

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
condition(du, t, integrator, N, M) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 31

Temp_rich = range(0, num_temps-1, length = num_temps)

progress = Progress((length(path)-1)*num_temps; desc="Progress running:")

all_Rrela_collect = Vector{Vector{Float64}}(); all_Crela_collect = Vector{Vector{Float64}}();
all_R_collect = Vector{Vector{Float64}}(); all_C_collect = Vector{Vector{Float64}}();
all_ii_collect = Vector{Vector{Float64}}(); all_ij_collect = Vector{Vector{Float64}}();
all_ii_sur_collect = Vector{Vector{Float64}}(); all_ij_sur_collect = Vector{Vector{Float64}}(); 
all_r_collect = Vector{Vector{Float64}}(); all_r_sur_collect = Vector{Vector{Float64}}()
Eff_results = zeros(Float64, num_temps, 65)
@time for j in 1: num_temps
    all_ℵii_H = Float64[]; all_ℵij_H = Union{Float64, Missing}[]; all_ℵij_d_H = Union{Float64, Missing}[];
    all_uℵij_H = Union{Float64, Missing}[]; all_lℵij_H = Union{Float64, Missing}[]; 
    all_ℵii_sur_H = Union{Float64, Missing}[]; all_ℵij_sur_H = Union{Float64, Missing}[]; all_ℵij_d_sur_H = Union{Float64, Missing}[];
    all_uℵij_sur_H = Union{Float64, Missing}[]; all_lℵij_sur_H = Union{Float64, Missing}[];
    all_r_H = Float64[]; all_r_sur_H = Float64[];
    all_leading_H = Float64[]; all_diag_H = Float64[];radi_H = Float64[]; diag_dominance_H = Float64[];
    all_u_H = Float64[]; all_m_H = Float64[]; 
    RO_H = Float64[]; ulO_H = Float64[]; Rul_H = Float64[]; 
    RO_sur_H = Union{Float64, Missing}[]; ulO_sur_H = Union{Float64, Missing}[];  Rul_sur_H = Union{Float64, Missing}[];
    all_Eu_H = Float64[]; all_Em_H = Float64[]; all_Eu_sur_H = Float64[]; all_Em_sur_H = Float64[];
    all_Tpu_H = Float64[]; all_Tpm_H = Float64[]; all_Tpu_sur_H = Float64[]; all_Tpm_sur_H = Float64[]; all_sumαij = Float64[];
    all_Rrela_H = Float64[]; all_Crela_H = Float64[]; all_R_H =  Float64[]; all_C_H = Float64[]

    for i in 1:length(path)
        # if i != 110 # need to remove this in p-1_new
        # if i != 485 && i != 211 && i != 678 # need to remove this in p0_new
            @load path[i] all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_ℵii_sur all_ℵij_sur all_ℵij_d_sur all_uℵij_sur all_lℵij_sur all_r all_r_sur all_u all_m RO ulO Rul RO_sur ulO_sur Rul_sur all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_diag radi diag_dominance all_Rrela all_Crela all_R all_C
            append!(all_ℵii_H, all_ℵii[j]); append!(all_ℵij_H, all_ℵij[j]); append!(all_ℵij_d_H,all_ℵij_d[j]); append!(all_uℵij_H,all_uℵij[j]); append!(all_lℵij_H,all_lℵij[j]);
            append!(all_ℵii_sur_H, all_ℵii_sur[j]); append!(all_ℵij_sur_H, all_ℵij_sur[j]); append!(all_ℵij_d_sur_H, all_ℵij_d_sur[j]); append!(all_uℵij_sur_H, all_uℵij_sur[j]); append!(all_lℵij_sur_H, all_lℵij_sur[j]);
            append!(all_r_H, all_r[j]); append!(all_r_sur_H, all_r_sur[j]);
            append!(all_u_H, all_u[j]); append!(all_m_H, all_m[j]); 
            append!(RO_H, RO[j]); append!(ulO_H, ulO[j]); append!(Rul_H, Rul[j]);
            append!(RO_sur_H, RO_sur[j]); append!(ulO_sur_H, ulO_sur[j]); append!(Rul_sur_H, Rul_sur[j]);
            append!(all_Eu_H, all_Eu[j]); append!(all_Em_H, all_Em[j]); append!(all_Eu_sur_H, all_Eu_sur[j]); append!(all_Em_sur_H, all_Em_sur[j]);
            append!(all_Tpu_H, all_Tpu[j]); append!(all_Tpm_H, all_Tpm[j]); append!(all_Tpu_sur_H, all_Tpu_sur[j]); append!(all_Tpm_sur_H, all_Tpm_sur[j]);
            push!(all_leading_H, all_leading[j]); append!(all_diag_H, all_diag[j]); append!(radi_H, radi[j]); push!(diag_dominance_H, diag_dominance[j]);
            append!(all_Rrela_H, all_Rrela[j]); append!(all_Crela_H, all_Crela[j]); append!(all_R_H, all_R[j]); append!(all_C_H, all_C[j])
        # end 
        next!(progress)
    end 

    push!(all_Rrela_collect, all_Rrela_H); push!(all_Crela_collect,all_Crela_H);
    push!(all_R_collect, all_R_H); push!(all_C_collect, all_C_H);
    push!(all_ii_collect, all_ℵii_H); push!(all_ij_collect, all_ℵij_H);
    push!(all_ii_sur_collect, all_ℵii_sur_H); push!(all_ij_sur_collect, all_ℵij_sur_H);
    push!(all_r_collect, all_r_H); push!(all_r_sur_collect, all_r_sur_H); 
    Eff_results[Int(j),:] = [mean(all_ℵii_H), std(all_ℵii_H)/sqrt(length(all_ℵii_H)), mean(skipmissing(all_ℵij_H)), std(skipmissing(all_ℵij_H))/sqrt(length(all_ℵij_H)), 
    mean(skipmissing(all_ℵij_d_H)), std(skipmissing(all_ℵij_d_H))/sqrt(length(all_ℵij_d_H)), mean(skipmissing(all_uℵij_H)), std(skipmissing(all_uℵij_H))/sqrt(length(all_uℵij_H)), mean(skipmissing(all_lℵij_H)), std(skipmissing(all_lℵij_H))/sqrt(length(all_lℵij_H)), 
    mean(all_ℵii_sur_H), std(all_ℵii_sur_H)/sqrt(length(all_ℵii_sur_H)), mean(skipmissing(all_ℵij_sur_H)), std(skipmissing(all_ℵij_sur_H))/sqrt(length(all_ℵij_sur_H)), 
    mean(skipmissing(all_ℵij_d_sur_H)), std(skipmissing(all_ℵij_d_sur_H))/sqrt(length(all_ℵij_d_sur_H)), mean(skipmissing(all_uℵij_sur_H)), std(skipmissing(all_uℵij_sur_H))/sqrt(length(all_uℵij_sur_H)), mean(skipmissing(all_lℵij_sur_H)), std(skipmissing(all_lℵij_sur_H))/sqrt(length(all_lℵij_sur_H)), 
    mean(all_r_H), std(all_r_H)/sqrt(length(all_r_H)), mean(all_r_sur_H), std(all_r_sur_H)/sqrt(length(all_r_sur_H)),
    mean(all_u_H), std(all_u_H)/sqrt(length(all_u_H)), mean(all_m_H), std(all_m_H)/sqrt(length(all_m_H)), 
    mean(skipmissing(RO_H)), std(skipmissing(RO_H))/sqrt(length(RO_H)),
    mean(skipmissing(ulO_H)), std(skipmissing(ulO_H))/sqrt(length(ulO_H)),
    mean(skipmissing(Rul_H)), std(skipmissing(Rul_H))/sqrt(length(Rul_H)),
    mean(skipmissing(RO_sur_H)), std(skipmissing(RO_sur_H))/sqrt(length(RO_sur_H)),
    mean(skipmissing(ulO_sur_H)), std(skipmissing(ulO_sur_H))/sqrt(length(ulO_sur_H)),
    mean(skipmissing(Rul_sur_H)), std(skipmissing(Rul_sur_H))/sqrt(length(Rul_sur_H)),
    mean(all_Eu_H), std(all_Eu_H)/sqrt(length(all_Eu_H)), mean(all_Em_H), std(all_Em_H)/sqrt(length(all_Em_H)),
    mean(all_Eu_sur_H), std(all_Eu_sur_H)/sqrt(length(all_Eu_sur_H)), mean(all_Em_sur_H), std(all_Em_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_Tpu_H), std(all_Tpu_H)/sqrt(length(all_Tpu_H)), mean(all_Tpm_H), std(all_Tpm_H)/sqrt(length(all_Tpm_H)),
    mean(all_Tpu_sur_H), std(all_Tpu_sur_H)/sqrt(length(all_Tpm_sur_H)), mean(all_Tpm_sur_H), std(all_Tpm_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_leading_H), std(all_leading_H)/sqrt(length(all_leading_H)), sum(real.(all_leading_H) .< 0)/length(path), 
    mean(diag_dominance_H), std(diag_dominance_H)/sqrt(length(diag_dominance_H)), 
    mean(all_diag_H), std(all_diag_H)/sqrt(length(all_diag_H)), mean(radi_H), std(radi_H)/sqrt(length(radi_H))]
end 
R"library(beepr); beep(sound = 4, expr = NULL)"
col_names_EF = ["αii", "αii_err", "αij", "αij_err", "αij_d", "αij_d_err", "αij_upper", "αij_upper_err","αij_lower", "αij_lower_err",
                "αii_sur", "αii_sur_err", "αij_sur", "αij_sur_err", "αij_d_sur", "αij_d_sur_err", "αij_upper_sur", "αij_upper_sur_err","αij_lower_sur", "αij_lower_sur_err",
                "r", "r_err", "r_sur", "r_sur_err", 
                "u", "u_err","m", "m_err", 
                "RO", "RO_err", "ulO", "ulO_err", "estα", "estα_err",
                "RO_sur", "RO_sur_err", "ulO_sur", "ulO_sur_err", "estα_sur", "estα_sur_err",
                "Eu", "Eu_err", "Em", "Em_err", "Eu_sur", "Eu_sur_err", "Em_sur", "Em_sur_err",
                "Tpu", "Tpu_err", "Tpm", "Tpm_err", "Tpu_sur", "Tpu_sur_err", "Tpm_sur", "Tpm_sur_err",
                "eigen", "eigen_err", "stability",
                "diag_dom", "diag_dom_err",
                "Jac_diag", "Jac_diag_err", "radius", "radius_err"];
Eff_results = DataFrame(Eff_results, col_names_EF);

##### Saving #######
# CSV.write("../results/Eff_results_p-1_new.csv", Eff_results, writeheader=false)
# @save "../results/Feas_CR_dist_p-1_new.jld2" all_Rrela_collect all_Crela_collect all_R_collect all_C_collect all_ii_collect all_ij_collect all_ii_sur_collect all_ij_sur_collect all_r_collect all_r_sur_collect

##### Loading #######
Eff_results = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)
@load "../results/Feas_CR_dist_p-1_new.jld2" all_Rrela_collect all_Crela_collect all_R_collect all_C_collect all_ii_collect all_ij_collect all_ii_sur_collect all_ij_sur_collect all_r_collect all_r_sur_collect

cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])
f = Figure(fontsize = 30, size = (1800, 600));
Label(f[:,0], "Minimal Trade-off", fontsize = 50, rotation = pi/2)
ax1 = Axis(f[1,1],# title = "Maximal ", # 
    xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
# ax2 = Axis(f[1,1], ylabel = "log(|αij/αii|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αii, color = ("#FA8328", 0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, Eff_results.αii .- Eff_results.αii_err, Eff_results.αii .+ Eff_results.αii_err, color = ("#FA8328", 0.2))
# lines!(ax1, Temp_rich, Eff_results.αii_sur, color = ("#FA8328", 1.0), linewidth = 5, label = "survivor αii")
# band!(ax1, Temp_rich, Eff_results.αii_sur .- Eff_results.αii_sur_err, Eff_results.αii_sur .+ Eff_results.αii_sur_err, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, Eff_results.αij, color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich,  Eff_results.αij .- Eff_results.αij_err, Eff_results.αij .+ Eff_results.αij_err, color = ("#015845", 0.2))
# lines!(ax1, Temp_rich,Eff_results.αij_sur, color = ("#015845", 1.0), linewidth = 5, label = "survivor αij")
# band!(ax1, Temp_rich, Eff_results.αij_sur .- Eff_results.αij_sur_err, Eff_results.αij_sur .+ Eff_results.αij_sur_err, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[1,1, TopLeft()], "(a)")
ax2 = Axis(f[1,2], limits = ((-22.0, 3.0), nothing), xlabel = "log(|αᵢᵢ|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
for i in 1: 31
    density!(ax2, log.(abs.(all_ii_collect[i])), color = (cs[i], 0.2), strokewidth = 3, strokecolor = (cs[i], 0.7))
    # lines!(ax, [median(log.(abs.(vcat(all_ℵii[i], all_ℵij[i])))), median(log.(abs.(vcat(all_ℵii[i], all_ℵij[i]))))], [0, 0.7], linestyle = :dash, color = (cs[i], 0.3), linewidth = 3)
end 
Label(f[1,2, TopLeft()], "(b)")
ax3 = Axis(f[1,3], limits = ((-22.0, 3.0), nothing), xlabel = "log(|αᵢⱼ|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
for i in 1: 31
    density!(ax3, log.(abs.(all_ij_collect[i])), color = (cs[i], 0.3), strokewidth = 3, strokecolor = (cs[i], 0.5))
    # lines!(ax, [median(log.(abs.(vcat(all_ℵii[i], all_ℵij[i])))), median(log.(abs.(vcat(all_ℵii[i], all_ℵij[i]))))], [0, 0.7], linestyle = :dash, color = (cs[i], 0.3), linewidth = 3)
end 
Colorbar(f[1,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[1,3, TopLeft()], "(c)")
f
save("../results/distα0.pdf", f) 

repeating_ii = [vcat([repeat([all_ii_collect[t][i]], 99) for i in 1:length(all_ii_collect[t])]...) for t in 1:num_temps]
cor_ij_ii = [cor(repeating_ii[t], all_ij_collect[t]) for t in 1:num_temps]

plot([mean(all_ij_collect[i]) for i in 1:num_temps])

plot(cor_ij_ii)

#######ONLY collecting CR
path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_C_collect = Vector{Vector{Float64}}(); all_R_collect = Vector{Vector{Float64}}()
all_Rrela_collect = Vector{Vector{Float64}}(); all_Crela_collect = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_R_H = Float64[]; all_C_H = Union{Float64, Missing}[];
    all_Rrela_H = Float64[]; all_Crela_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_R all_C all_Rrela all_Crela
        if sum(vcat(all_R...).> 1000) == 0
            append!(all_R_H,  all_R[j]); append!(all_C_H, all_C[j]);
            append!(all_Rrela_H,  all_Rrela[j]); append!(all_Crela_H, all_Crela[j])
        end 
        next!(progress)
    end 
    push!(all_R_collect, all_R_H); push!(all_C_collect, all_C_H);
    push!(all_Rrela_collect, all_Rrela_H); push!(all_Crela_collect, all_Crela_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

Temp_rich = range(0, num_temps-1, length = num_temps)
temp = collect(Temp_rich .+ 273.15)
temp_R = vcat([repeat([Temp_rich[t]], length(all_R_collect[t])) for t in 1:num_temps]...)
temp_C = vcat([repeat([Temp_rich[t]], length(all_C_collect[t])) for t in 1:num_temps]...)
allR = vcat(all_R_collect...)
allC = vcat(all_C_collect...)

# scatter(temp_R, allR)
meanR = [mean(all_R_collect[t]) for t in 1:num_temps]
varR = [var(all_R_collect[t]) for t in 1: num_temps]
R_err = [std(all_R_collect[t])/sqrt(length(all_R_collect[t])) for t in 1:num_temps]
meanC = [mean(all_C_collect[t]) for t in 1:num_temps]
varC = [var(all_C_collect[t])./mean(all_C_collect[t]) for t in 1: num_temps]
C_err = [std(all_C_collect[t])/sqrt(length(all_C_collect[t])) for t in 1:num_temps]


all_temp_R = vcat([repeat([temp[t]], length(all_Rrela_collect[t])) for t in 1:num_temps]...) .- 273.15
all_relaR = vcat(all_Rrela_collect...)
[mean(all_Rrela_collect[t]) for t in 1:num_temps]
f = Figure(fontsize = 35, size = (1200, 1200));
ax1 = Axis(f[1,1], xlabel = "Temperature", ylabel = "Resource distribution", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
# scatter!(ax1,all_temp_R .- 273.15, vcat(all_Rrela_collect...), color = ("#285C93"), label = "", alpha = 0.2)
boxplot!(ax1, all_temp_R, all_relaR)
f
# save("../results/Feas_R_dist_p-1.png", f) 
# println([mean(all_Rrela_collect[t]) for t in 1: num_temps] ,"\n")

all_temp_C = vcat([repeat([temp[t]], length(all_Crela_collect[t])) for t in 1:num_temps]...) .- 273.15
all_relaC = vcat(all_Crela_collect...)
f = Figure(fontsize = 35, size = (1200, 1200));
ax2 = Axis(f[1,1], xlabel = "Temperature", ylabel = "Consumer distribution", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
# scatter!(ax2,all_temp_C .- 273.15, vcat(all_Crela_collect...), color = ("#285C93"), label = "", alpha = 0.7)
boxplot!(ax2, all_temp_C, all_relaC)
f
# println([mean(all_Crela_collect[t]) for t in 1: num_temps] ,"\n")
# save("../results/Feas_C_dist_p-1.png", f) 


# var_αii = [var(all_ii_collect[t])/abs.(mean(all_ii_collect[t])) for t in 1:num_temps]
# var_αij = [var(all_ij_collect[t])/abs.(mean(all_ij_collect[t])) for t in 1:num_temps]
# plot(var_αii)

allαii = vcat(all_ii_collect...)
allαij = vcat(all_ij_collect...)

allα = vcat(allαii, allαij)

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "α", ylabel = "frequency (αii)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "frequency (αij)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hist!(ax1, allαii, color =( "#FA8328", 0.7), bins = 100)
hist!(ax2, allαij, color = ("#069F66", 0.7), bins = 100)
linkxaxes!(ax1,ax2)
p1 = [PolyElement(color = ("#FA8328", 0.7), strokecolor = :transparent)]
p2 = [PolyElement(color = ("#069F66", 0.7), strokecolor = :transparent)]
Legend(f[1,1], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
Label(f[1,1, TopLeft()], "(d)")
f
save("../results/a-1_T_hist.png", f) 


temp = collect(Temp_rich .+273.15)

temp_all_ii = vcat([repeat([Temp_rich[t]], length(all_ii_collect[t])) for t in 1:num_temps]...)
temp_all_ij = vcat([repeat([Temp_rich[t]], length(all_ij_collect[t])) for t in 1:num_temps]...)
temp_r = vcat([repeat([Temp_rich[t]], length(all_r_collect[t])) for t in 1:num_temps]...)

allαii = vcat(all_ii_collect...)
allαij = vcat(all_ij_collect...)
allr = vcat(all_r_collect...)

allα = vcat(allαii, allαij)
100* sum(allα .> 0)/length(allα)

scatter(temp_all_ii, allαii)
scatter(temp_all_ij, allαij)
scatter(temp_r, allr)




############## all α ################
path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_ij_collect = Vector{Vector{Float64}}() ; all_sym_ijii_collect = Vector{Vector{Float64}}(); all_sym_collect = Vector{Vector{Float64}}()
# all_αijii_collect = Vector{Vector{Float64}}()
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_ij_H = Float64[]; all_sym_ii_H = Float64[]; all_sym_H = Float64[]
    # all_αijii_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij #all_ℵii_sur
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        ij = [A[i, j]/A[i, i] for i in 1:N for j in 1:N if j != i]
        sym_ij = std([norm(A[i, :]) for i in 1:N])
        sym_ijii = [(A[i, j]/A[i, i] - A[j, i]/A[j, j])^2 / ((A[i, j]/A[i, i])^2 + (A[j, i]/A[j, j])^2) for i in 1:N for j in 1:N if j > i]
        αijii = var([A[i, j]/A[i, i] for i in 1:N for j in 1:N if j != i])
        append!(all_ij_H, ij); append!(all_sym_ii_H,  sym_ijii); append!(all_sym_H, sym_ij)
        next!(progress)
    end 
    push!(all_ij_collect, all_ij_H); push!(all_sym_ijii_collect, all_sym_ii_H); push!(all_sym_collect, all_sym_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

# @save "../results/Asymmetry0.jld2" all_sym_collect all_sym_ijii_collect # all_αijii_collect
# @load "../results/Asymmetry-1.jld2" all_sym_collect all_sym_ijii_collect all_αijii_collect

sym = [mean(all_sym_collect[t]) for t in 1: num_temps]
sym_err = [std(all_sym_collect[t])/sqrt(length(all_sym_collect[t])) for t in 1: num_temps]
symii = [mean(all_sym_ijii_collect[t]) for t in 1: num_temps]
symii_err = [std(all_sym_ijii_collect[t])/sqrt(length(all_sym_ijii_collect[t])) for t in 1: num_temps]
plot(Temp_rich, sym)
# f = Figure(fontsize = 35, size = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Variation log(αij)", ylabel = "Stability", xlabelsize = 50, ylabelsize = 50)
# scatter!(ax, var_ij, Eff_results.stability, color = "#015845", markersize = 12, alpha = 0.8)
# for (x, y, e) in zip(var_ij, Eff_results.stability, var_ij_err)
#     # Vertical line
#     lines!(ax, [x - e, x + e], [y, y], color = ("#015845", 0.4), linewidth = 1)
#     # Horizontal caps
#     cap_length = 0.001*mean(Eff_results.stability)  # Length of horizontal caps
#     # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
#     # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
#     lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = ("#015845", 0.4), linewidth = 1)
#     lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = ("#015845", 0.4), linewidth = 1)
# end
# Label(f[1,1, TopLeft()], "(b)")
# f
# save("../results/varij-1.pdf", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Asymmetry of α", xlabelsize = 50, ylabelsize = 50)
# ax2 = Axis(f[1,1], ylabel = "log(|α|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich[1:30], symii[1:30], color = ("#285C93", 0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, symii - symii_err, symii+ symii_err, color = ("#285C93", 0.4))
# # axislegend(position = :rb)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/sym0.pdf", f) 

############ rich vs. stability ############
path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_rich_collect_0 = Vector{Vector{Float64}}(); var_u_collect_0 = Vector{Vector{Float64}}(); var_ϵ_collect_0 = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_rich_H = Float64[]; var_u_H = Float64[]; var_ϵ_H = Float64[]
    for i in 1:length(path)
        @load path[i] all_ℵii_sur all_u all_m
        var_u = std(log.(abs.(all_u[j])))
        ϵ = (all_u[j] .* (1 - 0.3) .- all_m[j]) ./ (all_u[j])
        append!(var_u_H, var_u); append!(all_rich_H, length(all_ℵii_sur[j]))
        append!(var_ϵ_H, std(ϵ))
        next!(progress)
    end 
    push!(all_rich_collect_0, all_rich_H); push!(var_u_collect_0, var_u_H); push!(var_ϵ_collect_0, var_ϵ_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

stability_0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 59]
stability_1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 59]
mean_rich_0 = [mean(all_rich_collect_0[t]) for t in 1:num_temps]
rich_err_0 = [std(all_rich_collect_0[t])/sqrt(length(all_rich_collect_0[t])) for t in 1: num_temps]
mean_u_0 = [mean(var_u_collect_0[t]) for t in 1:num_temps]
u_err_0 = [std(var_u_collect_0[t])/sqrt(length(var_u_collect_0[t])) for t in 1: num_temps]
mean_rich_1 = [mean(all_rich_collect_1[t]) for t in 1:num_temps]
rich_err_1 = [std(all_rich_collect_1[t])/sqrt(length(all_rich_collect_1[t])) for t in 1: num_temps]
mean_u_1 = [mean(var_u_collect_1[t]) for t in 1:num_temps]
u_err_1 = [std(var_u_collect_1[t])/sqrt(length(var_u_collect_1[t])) for t in 1: num_temps]

mean_ϵ_0 = [mean(log.(var_ϵ_collect_0[t])) for t in 1:num_temps]
mean_ϵ_1 = [mean(log.(var_ϵ_collect_1[t])) for t in 1:num_temps]

stability_all = vcat(stability_0, stability_1)
mean_rich = vcat(mean_rich_0, mean_rich_1)
rich_err = vcat(rich_err_0, rich_err_1)
u_all = vcat(mean_u_0, mean_u_1)
ϵ_all = vcat(mean_ϵ_0, mean_ϵ_1)
# u_all_err = vcat(u_err_0, u_err_1)
cs = ColorScheme(range(colorant"#F5D44B",colorant"#20090E",  length = length(ϵ_all)))
f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Richness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
sc = scatter!(ax, mean_rich, stability_all, color = ϵ_all, colormap = cs, markersize = 15, alpha = 0.8)
for (x, y, e) in zip(mean_rich, stability_all, rich_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#F5D44B", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(stability_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = ("#F5D44B", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = ("#F5D44B", 0.4), linewidth = 1)
end
subgrid = GridLayout(f[1,2], tellheight = false)
Colorbar(f[1,2], colorrange = [minimum(ϵ_all), maximum(ϵ_all)], colormap = cs, label = "σ(log(u))")
f



f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Richness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, mean_rich, stability_all, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(mean_rich, stability_all, rich_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(stability_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = ("#376298", 0.4), linewidth = 1)
end
# Label(f[1,1, TopLeft()], "(a)")
f

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Richness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, mean_rich_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(mean_rich_0, stability_0, rich_err_0)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, mean_rich_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(mean_rich_1, stability_1, rich_err_1)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rb)
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/rich_sta.pdf", f) 



######## 
diag_dom_0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:,60]
diag_dom_0_err = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:,61]
diag_dom_1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:,60]
diag_dom_1_err = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:,61]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "p(Diagonal Dominance)", ylabel = "Richness", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, diag_dom_0, mean_rich_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e_r, e_s) in zip(diag_dom_0, mean_rich_0, rich_err_0, diag_dom_0_err)
    # Vertical line
    lines!(ax, [x - e_s, x + e_s], [y, y], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x , x], [y-e_r, y+e_r], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.0001 * mean(diag_dom_0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e_s, x - e_s], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e_s, x + e_s], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y - e_r, y - e_r], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y + e_r, y + e_r], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, diag_dom_1, mean_rich_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e_r, e_s) in zip(diag_dom_1, mean_rich_1, rich_err_1, diag_dom_1_err)
    # Vertical line
    lines!(ax, [x - e_s, x + e_s], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x , x], [y-e_r, y+e_r], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.0001 * mean(diag_dom_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e_s, x - e_s], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e_s, x + e_s], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y - e_r, y - e_r], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y + e_r, y + e_r], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rt)
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/diag_dom_rich.pdf", f) 


################# evenness ##################
function pielou_evenness(C_rela)
    H = -sum(C_rela .* log.(C_rela))
    H_max = log(length(C_rela))
    return H / H_max
end

path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_even_collect_0 = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_even_H = Float64[];
    # all_αijii_H = Float64[]
    for i in 1:length(path)
        @load path[i] all_Crela
        even = pielou_evenness(all_Crela[j])
        append!(all_even_H, even)
        # append!(all_αijii_H,  αijii);
        next!(progress)
    end 
    push!(all_even_collect_0, all_even_H)
    # push!(all_αijii_collect, all_αijii_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

even_mean_0 = [mean(all_even_collect_0[t]) for t in 1: num_temps]
even_err_0 = [std(all_even_collect_0[t])/sqrt(length(all_even_collect_0[t])) for t in 1: num_temps]
even_mean_1 = [mean(all_even_collect_1[t]) for t in 1: num_temps]
even_err_1 = [std(all_even_collect_1[t])/sqrt(length(all_even_collect_1[t])) for t in 1: num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Evenness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, even_mean_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(even_mean_0, stability_0, even_err_0)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, even_mean_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(even_mean_1, stability_1, even_err_1)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rb)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/sta_varC.pdf", f) 


even_all = vcat(even_mean_0, even_mean_1)
var_u_all = vcat(mean_u_0, mean_u_1)

even_mean_0 = [mean(all_even_collect_0[t]) for t in 1: num_temps]
even_err_0 = [std(all_even_collect_0[t])/sqrt(length(all_even_collect_0[t])) for t in 1: num_temps]
even_mean_1 = [mean(all_even_collect_1[t]) for t in 1: num_temps]
even_err_1 = [std(all_even_collect_1[t])/sqrt(length(all_even_collect_1[t])) for t in 1: num_temps]

mean_ϵ_0 = [mean(log.(var_ϵ_collect_0[t])) for t in 1:num_temps]
ϵ_err_0 = [std(log.(var_ϵ_collect_0[t]))/sqrt(length(var_ϵ_collect_0[t])) for t in 1: num_temps]
mean_ϵ_1 = [mean(log.(var_ϵ_collect_1[t])) for t in 1:num_temps]
ϵ_err_1 = [std(log.(var_ϵ_collect_1[t]))/sqrt(length(var_ϵ_collect_1[t])) for t in 1: num_temps]

ϵ_all = vcat(mean_ϵ_0, mean_ϵ_1)
ϵ_err = vcat(ϵ_err_0, ϵ_err_1)

rich_all = vcat(mean_rich_0, mean_rich_1)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(ϵ)", ylabel = "Evenness", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, mean_ϵ_0, even_mean_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e_r, e_s) in zip(mean_ϵ_0, even_mean_0, ϵ_err_0, even_err_0)
    # Vertical line
    lines!(ax, [x - e_s, x + e_s], [y, y], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x , x], [y-e_r, y+e_r], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.0001 * mean(even_mean_0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e_s, x - e_s], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e_s, x + e_s], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y - e_r, y - e_r], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y + e_r, y + e_r], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, mean_ϵ_1, even_mean_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e_r, e_s) in zip(mean_ϵ_1, even_mean_1, ϵ_err_1, even_err_1)
    # Vertical line
    lines!(ax, [x - e_s, x + e_s], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x , x], [y-e_r, y+e_r], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.0001 * mean(even_mean_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e_s, x - e_s], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e_s, x + e_s], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y - e_r, y - e_r], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_0, x + cap_length_0], [y + e_r, y + e_r], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rt)
# Label(f[1,1, TopLeft()], "(a)")
f


############### σ vs. stability 
stability_0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 59]
stability_1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 59]

path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_ij_collect_1 = Vector{Vector{Float64}}(); all_ii_collect_1 = Vector{Vector{Float64}}();
all_ij_collect_org_1 = Vector{Vector{Float64}}(); all_ii_collect_org_1 = Vector{Vector{Float64}}();
var_u_collect_1 = Vector{Vector{Float64}}() ; var_m_collect_1 = Vector{Vector{Float64}}()
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_ij_H = Float64[]; all_ii_H = Float64[]; 
    all_ij_org_H = Float64[]; all_ii_org_H = Float64[]; 
    var_u_H = Float64[]; var_m_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij all_u all_m #all_ℵii_sur
        # A = zeros(Float64, N, N)
        # A[ind_off] = all_ℵij[j]
        # A[diagind(A)] = all_ℵii[j]
        # ij = std(log.(abs.([A[i, j]/A[i, i] for i in 1:N for j in 1:N if j != i])))
        ij = std(log.(abs.(all_ℵij[j])))
        ii = std(log.(abs.(all_ℵii[j])))
        var_u = std(log.(abs.(all_u[j])))
        var_m = std(log.(abs.(all_m[j])))
        append!(all_ij_H, ij); append!(all_ii_H, ii); 
        append!(all_ij_org_H, all_ℵij[j]); append!(all_ii_org_H, all_ℵii[j]); 
        append!(var_u_H, var_u); append!(var_m_H, var_m); 
        next!(progress)
    end 
    push!(all_ij_collect_1, all_ij_H); push!(all_ii_collect_1, all_ii_H); 
    push!(all_ij_collect_org_1, all_ij_org_H); push!(all_ii_collect_org_1, all_ii_org_H); 
    push!(var_u_collect_1, var_u_H); push!( var_m_collect_1, var_m_H); 
end 

all_ii = vcat(all_ii_collect_0...,all_ii_collect_1...)
all_ij = vcat(all_ij_collect_0...,all_ij_collect_1...)
all_u = vcat(var_u_collect_0...,var_u_collect_1...)
all_m = vcat(var_m_collect_0...,var_m_collect_1...)

function ellipse_points(center, axis_lengths, angle_ellipse; n_points=100)
    t = range(0, 2π, length=n_points)
    ellipse = [axis_lengths[1] .* cos.(t), axis_lengths[2] .* sin.(t)]  
    rotation_matrix = [cos(angle_ellipse) -sin(angle_ellipse); sin(angle_ellipse) cos(angle_ellipse)]
    rotated_ellipse = rotation_matrix * ellipse
    [rotated_ellipse[i] .+ center[i] for i in 1:2] 
end
confidence_interval_factor = 2.4477  # 95% confidence for 2D ellipse
f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(log(u))", ylabel = "σ(log(α))", xlabelsize = 50, ylabelsize = 50)
x = all_u; y = all_ii
mean_x, mean_y = mean(x), mean(y)
eigen_vals, eigen_vecs = eigen(cov(hcat(all_u, all_ii)))
axis_lengths = sqrt.(eigen_vals) * confidence_interval_factor
angle_ellipse = atan(eigen_vecs[2, 1], eigen_vecs[1, 1])
ellipse_pts_1 = ellipse_points([mean(all_u), mean(all_ii)], axis_lengths, angle_ellipse)
scatter!(ax, x, y, color = "#FA8328", markersize = 15, alpha = 0.1, label = "αᵢᵢ")
lines!(ax, ellipse_pts_1[1], ellipse_pts_1[2], color="#0758AE", linewidth=5)
x = all_u; y = all_ij
mean_x, mean_y = mean(x), mean(y)
eigen_vals, eigen_vecs = eigen(cov(hcat(all_u, all_ij)))
axis_lengths = sqrt.(eigen_vals) * confidence_interval_factor
angle_ellipse = atan(eigen_vecs[2, 1], eigen_vecs[1, 1])
ellipse_pts_2 = ellipse_points([mean(all_u), mean(all_ij)], axis_lengths, angle_ellipse)
scatter!(ax, x, y, color = "#015845", markersize = 15, alpha = 0.1, label = "αᵢⱼ")
lines!(ax, ellipse_pts_2[1], ellipse_pts_2[2], color="#0758AE", linewidth=5)
axislegend(position = :lt)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/σu_σα.pdf", f) 




# For αᵢᵢ
eigen_vals, eigen_vecs = eigen(cov(hcat(all_u, all_ii)))
axis_lengths = sqrt.(eigen_vals) * confidence_interval_factor
angle_ellipse = atan(eigen_vecs[2, 1], eigen_vecs[1, 1])
ellipse_pts_1 = ellipse_points([mean(all_u), mean(all_ii)], axis_lengths, angle_ellipse)

# For αᵢⱼ
eigen_vals, eigen_vecs = eigen(cov(hcat(all_u, all_ij)))
axis_lengths = sqrt.(eigen_vals) * confidence_interval_factor
angle_ellipse = atan(eigen_vecs[2, 1], eigen_vecs[1, 1])
ellipse_pts_2 = ellipse_points([mean(all_u), mean(all_ij)], axis_lengths, angle_ellipse)




# lm_model = lm(@formula(ii ~ u), df)
# coeffs = coef(lm_model)
# u_range = range(minimum(all_u), stop=maximum(all_u), length=100)
# fitted_line = coeffs[1] .+ coeffs[2] .* u_range
# plot!(u_range, fitted_line, label="Fitted Line", color=:red)

var_ij_0 = [mean(all_ij_collect_0[t]) for t in 1: num_temps]
var_err_0 = [std(all_ij_collect_0[t])/sqrt(length(all_ij_collect_0[t])) for t in 1: num_temps]
var_ii_0 = [mean(all_ii_collect_0[t]) for t in 1: num_temps]
var_ii_err_0 = [std(all_ii_collect_0[t])/sqrt(length(all_ii_collect_0[t])) for t in 1: num_temps]
var_u_0 = [mean(var_u_collect_0[t]) for t in 1: num_temps]
var_u_err_0 = [std(var_u_collect_0[t])/sqrt(length(var_u_collect_0[t])) for t in 1: num_temps]
var_m_0 = [mean(var_m_collect_0[t]) for t in 1: num_temps]
var_m_err_0 = [std(var_m_collect_0[t])/sqrt(length(var_m_collect_0[t])) for t in 1: num_temps]

var_ij_1 = [mean(all_ij_collect_1[t]) for t in 1: num_temps]
var_err_1 = [std(all_ij_collect_1[t])/sqrt(length(all_ij_collect_1[t])) for t in 1: num_temps]
var_ii_1 = [mean(all_ii_collect_1[t]) for t in 1: num_temps]
var_ii_err_1 = [std(all_ii_collect_1[t])/sqrt(length(all_ii_collect_1[t])) for t in 1: num_temps]
var_u_1 = [mean(var_u_collect_1[t]) for t in 1: num_temps]
var_u_err_1 = [std(var_u_collect_1[t])/sqrt(length(var_u_collect_1[t])) for t in 1: num_temps]
var_m_1 = [mean(var_m_collect_1[t]) for t in 1: num_temps]
var_m_err_1 = [std(var_m_collect_1[t])/sqrt(length(var_m_collect_1[t])) for t in 1: num_temps]

var_ij = vcat(var_ij_0, var_ij_1); var_ii = vcat(var_ii_0, var_ii_1);
stability_all = vcat(stability_0, stability_1); 
var_ij_err = vcat(var_err_0, var_err_1); var_ii_err = vcat(var_ii_err_0, var_ii_err_1); 

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(log(α))", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, var_ii, stability_all, color = "#FA8328", markersize = 15, alpha = 0.8, label = "αᵢᵢ")
for (x, y, e) in zip(var_ii, stability_all, var_ii_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#FA8328", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#FA8328", 0.4), linewidth = 1)
end
scatter!(ax, var_ij, stability_all, color = "#015845", markersize = 15, alpha = 0.8, label = "αᵢⱼ")
for (x, y, e) in zip(var_ij, stability_all, var_ij_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#015845", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#015845", 0.4), linewidth = 1)
end
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/sta_σα.pdf", f) 

over_zero_0 = [sum(all_ii_collect_org_0[t].> 1.0e-7)/length(all_ii_collect_org_0[t]) for t in 1:num_temps]
over_zero_1 = [sum(all_ii_collect_org_1[t].> 1.0e-7)/length(all_ii_collect_org_1[t]) for t in 1:num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "p(αii > 0)", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, over_zero_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
scatter!(ax, over_zero_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(c)")
f
save("../results/sta_αii.pdf", f) 

reactivity_0 = [sum(all_leadH_0[t] .>0)./length(path_0) for t in 1:num_temps]
reactivity_1 = [sum(all_leadH_1[t] .>0)./length(path_1) for t in 1:num_temps]
reactivity_all = vcat(reactivity_0, reactivity_1)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(log(α))", ylabel = "p(Reactivity)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, var_ii, reactivity_all, color = "#FA8328", markersize = 15, alpha = 0.8, label = "αᵢᵢ")
for (x, y, e) in zip(var_ii, reactivity_all, var_ii_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#FA8328", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.0001 * mean(reactivity_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#FA8328", 0.4), linewidth = 1)
end
scatter!(ax, var_ij, reactivity_all, color = "#015845", markersize = 15, alpha = 0.8, label = "αᵢⱼ")
for (x, y, e) in zip(var_ij, reactivity_all, var_ij_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#015845", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.0001 * mean(reactivity_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#015845", 0.4), linewidth = 1)
end
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(d)")
f
save("../results/react_σα.pdf", f) 


var_u = vcat(var_u_0, var_u_1)
var_u_err = vcat(var_u_err_0, var_u_err_1)
var_m = vcat(var_m_0, var_m_1)
var_m_err = vcat(var_m_err_0, var_m_err_1)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(log(u,m))", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, var_u, stability_all, color = "#015845", markersize = 15, alpha = 0.8, label = "u")
for (x, y, e) in zip(var_u, stability_all, var_u_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#015845", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#015845", 0.4), linewidth = 1)
end
scatter!(ax, var_m, stability_all, color = "#FA8328", markersize = 15, alpha = 0.8, label = "m")
for (x, y, e) in zip(var_m, stability_all, var_m_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#FA8328", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#FA8328", 0.4), linewidth = 1)
end
axislegend(position = :rt)
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/sta_σu&m.pdf", f) 

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(log(u))", ylabel = "σ(log(α))", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, var_u, var_ii, color = "#015845", markersize = 15, alpha = 0.8, label = "αᵢᵢ")
for (x, y, e_u, e_ii) in zip(var_u, var_ii, var_u_err, var_ii_err)
    # Vertical line
    lines!(ax, [x - e_u, x + e_u], [y, y], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x, x], [y - e_ii, y + e_ii], color = ("#015845", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(var_u)  # Length of horizontal caps
    lines!(ax, [x - e_u, x - e_u], [y - cap_length, y + cap_length], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x + e_u, x + e_u], [y - cap_length, y + cap_length], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length, x + cap_length], [y - e_ii, y - e_ii], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length, x + cap_length], [y + e_ii, y + e_ii], color = ("#015845", 0.4), linewidth = 1)
end
scatter!(ax, var_u, var_ij, color = "#FA8328", markersize = 15, alpha = 0.8, label = "αᵢⱼ")
for (x, y, e_u, e_ij) in zip(var_u, var_ij, var_u_err, var_ij_err)
    # Vertical line
    lines!(ax, [x - e_u, x + e_u], [y, y], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x, x], [y - e_ij, y + e_ij], color = ("#015845", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(var_u)  # Length of horizontal caps
    lines!(ax, [x - e_u, x - e_u], [y - cap_length, y + cap_length], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x + e_u, x + e_u], [y - cap_length, y + cap_length], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length, x + cap_length], [y - e_ij, y - e_ij], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length, x + cap_length], [y + e_ij, y + e_ij], color = ("#FA8328", 0.4), linewidth = 1)
end
axislegend(position = :lt)
# Label(f[1,1, TopLeft()], "(a)")
f

########### upper lower #############

path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_upper_collect = Vector{Vector{Float64}}() ; all_lower_collect = Vector{Vector{Float64}}(); all_dm_collect = Vector{Vector{Float64}}()
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_upper_H = Float64[]; all_lower_H = Float64[]; all_dm_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        A_off = [sum(abs.(A[i, j]) for j in 1:N if j != i) for i in 1:N ]
        diag_dom = sum(abs.(all_ℵii[j]) - A_off .> 0)/N
        upper = [A[i, j]/A[i, i] for i in 1:N for j in 1:N if j > i]
        lower = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j > i]
        append!(all_upper_H, upper); append!(all_lower_H,  lower); append!(all_dm_H,  diag_dom)
        next!(progress)
    end 
    push!(all_upper_collect, all_upper_H); push!(all_lower_collect, all_lower_H); push!(all_dm_collect, all_dm_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"


plot(log.(abs.(all_upper_collect[1])), log.(abs.(all_lower_collect[1])))

### only upper lower
path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_upper_collect = Vector{Vector{Float64}}() ; all_lower_collect = Vector{Vector{Float64}}()
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_upper_H = Float64[]; all_lower_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        A_off = [sum(abs.(A[i, j]) for j in 1:N if j != i) for i in 1:N ]
        diag_dom = sum(abs.(all_ℵii[j]) - A_off .> 0)/N
        upper = [A[i, j] for i in 1:N for j in 1:N if j > i]
        lower = [A[j, i] for i in 1:N for j in 1:N if j > i]
        append!(all_upper_H, upper); append!(all_lower_H,  lower)
        next!(progress)
    end 
    push!(all_upper_collect, all_upper_H); push!(all_lower_collect, all_lower_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

cor(vcat(all_upper_collect...), vcat(all_lower_collect...))


########### diagonal dominance ###########
diag_dom_0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:,60]
diag_dom_0_err = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:,61]
diag_dom_1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:,60]
diag_dom_1_err = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:,61]

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "p(Diagnoal Dominance)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "p(Diagnoal Dominance)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich,diag_dom_0, color = ("#376298", 0.8), linewidth = 5, label = "ρ = 0")
band!(ax1, Temp_rich, diag_dom_0 .- diag_dom_0_err, diag_dom_0 .+ diag_dom_0_err, color = ("#376298", 0.2))
lines!(ax2, Temp_rich, diag_dom_1, color = ("#9A2B1A", 0.8), linewidth = 5, label = "ρ = -1")
band!(ax2, Temp_rich, diag_dom_1 .- diag_dom_1_err, diag_dom_1 .+ diag_dom_1_err, color = ("#9A2B1A", 0.2))
l1 = [LineElement(color = ("#376298",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#9A2B1A", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, ["ρ = 0", "ρ = -1"], halign = :left, valign = :top)
# Label(f[1,1, TopLeft()], "(a)")
f

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "p(Diagnoal Dominance)", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, diag_dom_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(diag_dom_0, stability_0, diag_dom_0_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, diag_dom_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(diag_dom_1, stability_1, diag_dom_1_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :lb)
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/sta_diag_dom.pdf", f) 

##########
path_1 = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path_0)*num_temps; desc="Progress running:")
num_temps = 31
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
all_leading_collect1 = Vector{Vector{ComplexF64}}(); α_ij_collect1 = Vector{Vector{Float64}}(); α_ii_collect1 = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_leading_H = ComplexF64[]; α_ij = Float64[]; α_ii = Float64[]
    for i in 1:length(path_1)
        @load path_1[i] all_ℵii all_ℵij all_r_sur all_ℵii_sur all_C
        sur = findall(x -> x in all_ℵii_sur[j], all_ℵii[j])
        N_s = length(sur)
        C = all_C[j]; r = all_r_sur[j]
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        A_sur = A[sur, sur]
        LV_Jac = [A_sur[i, j]*C[i] for i in 1:N_s, j in 1:N_s]
        # LV_Jac[diagind(LV_Jac)] .= [r[i] + A_sur[i, i]*C[i] + sum(A_sur[i, j]*C[j] for j in 1:N_s) for i in 1:N_s]
        jac_eigen = eigen(LV_Jac).values
        leading = jac_eigen[argmax(real.(jac_eigen))]
        append!(all_leading_H,  leading); append!(α_ij, all_ℵij[j]); append!(α_ii, all_ℵii[j])
        next!(progress)
    end 
    push!(all_leading_collect1, all_leading_H); push!(α_ij_collect1, α_ij); push!(α_ii_collect1, α_ii)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"
sta0 = [sum(real.(all_leading_collect0[t]) .< 0)/length(path_0) for t in 1: num_temps]
sta1 = [sum(real.(all_leading_collect1[t]) .< 0)/length(path_1) for t in 1: num_temps]

mean_rich_0 = [mean(all_rich_collect_0[t]) for t in 1:num_temps]
rich_err_0 = [std(all_rich_collect_0[t])/sqrt(length(all_rich_collect_0[t])) for t in 1: num_temps]
mean_rich_1 = [mean(all_rich_collect_1[t]) for t in 1:num_temps]
rich_err_1 = [std(all_rich_collect_1[t])/sqrt(length(all_rich_collect_1[t])) for t in 1: num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Richness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, mean_rich_0, sta0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(mean_rich_0, sta0, rich_err_0)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(sta0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, mean_rich_1, sta1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(mean_rich_1, sta1, rich_err_1)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(sta1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rb)
# Label(f[1,1, TopLeft()], "(a)")
f


############## stability ################
path = glob("Eff_iters*", "../data/Eff_p0_new/")

idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
progress = Progress(length(path)*num_temps; desc="Progress running:")
all_circ_0 = Vector{Vector{ComplexF64}}(); all_ij_collect_0 = Vector{Vector{Float64}}()
for j in 1: num_temps
    # if (j-1) % 5 == 0
        circ_leading_H = ComplexF64[]; all_ij_H = Float64[]
        for i in 1:length(path)
            @load path[i] all_leading all_ℵij all_ℵii
            A = zeros(Float64, N, N)
            A[ind_off] = all_ℵij[j]
            A[diagind(A)] = all_ℵii[j]
            ij = std([A[i, j]/A[i, i] for i in 1:N for j in 1:N if j != i])
            next!(progress)
            push!(circ_leading_H, all_leading[j]); push!(all_ij_H, ij)
        end 
    # end 
    push!(all_circ_0, circ_leading_H); push!(all_ij_collect_0, all_ij_H)
end 

all_var = vcat(all_ij_collect_0..., all_ij_collect_1...) 
all_real = real.(vcat(all_circ_0..., all_circ_1...))
# var_ij_0 = [mean(all_ij_collect_0[t]) for t in 1: num_temps]
# var_err_0 = [std(all_ij_collect_0[t])/sqrt(length(all_ij_collect_0[t])) for t in 1: num_temps]
# var_ij_1 = [mean(all_ij_collect_1[t]) for t in 1: num_temps]
# var_err_1 = [std(all_ij_collect_1[t])/sqrt(length(all_ij_collect_1[t])) for t in 1: num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "σ(log(α))", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, var_ij_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(var_ij_0, stability_0, var_err_0)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_0)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, var_ij_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(var_ij_1, stability_1, var_err_1)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rt)
# Label(f[1,1, TopLeft()], "(a)")
f


temp = hcat([repeat([Temp_rich[t]], length(all_circ[t])) for t in 1:num_temps]...)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Real λᴶₘ ", xlabelsize = 50, ylabelsize = 50)
for t in 1:num_temps
    scatter!(ax, temp[:,t], real.(all_circ[t]), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
end
lines!(ax, [0, 30], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/leading_ρ-1.pdf", f) 

path = glob("Eff_iters*", "../data/Eff_p0_ri/")

all_sta = Float64[]; all_leading_collect = Vector{Vector{ComplexF64}}()
for j in 1: num_temps
        circ_leading_H = ComplexF64[];
        for i in 1:length(path)
            @load path[i] all_leading
            push!(circ_leading_H, all_leading[j])
        end 
    push!(all_sta, sum(real.(circ_leading_H) .< 0)/length(path)); push!(all_leading_collect, circ_leading_H)
end 

########### Reactivity ############# 
path_0 = glob("Eff_iters*", "../data/Eff_p0_ri/")
all_leadH_0 = Vector{Vector{Float64}}(); all_circ_0 = Vector{Vector{ComplexF64}}()
for j in 1: num_temps
    # if (j-1) % 5 == 0
    all_leadH_H = Float64[]; circ_leading_H = ComplexF64[]
        for i in 1:length(path_0)
            @load path_0[i] all_H_leading all_leading
            push!(all_leadH_H, all_H_leading[j]); push!(circ_leading_H, all_leading[j])
        end 
    # end 
    push!(all_leadH_0, all_leadH_H); push!(all_circ_0, circ_leading_H)
end 
plot([sum(all_leadH_0[t] .>0)./length(path) for t in 1:num_temps],[sum(real.(all_circ_0[t]) .<0)./length(path) for t in 1:num_temps])

count(x -> x[1] < 0 && x[2] > 0, zip(real.(vcat(all_circ_0...)), vcat(all_leadH_0...)))/ length(circ_vec)

temp_0 = hcat([repeat([Temp_rich[t]], length(all_leadH_0[t])) for t in 1:num_temps]...)
temp_1 = hcat([repeat([Temp_rich[t]], length(all_leadH_1[t])) for t in 1:num_temps]...)

f = Figure(fontsize = 30, size = (1800, 600));
Label(f[:,0], "Minimal Trade-off", fontsize = 50, rotation = pi/2)
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Real λᴶₘ ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax1_2 = Axis(f[1,1], ylabel = "p(Stability)", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
for t in 1: num_temps
    scatter!(ax1, temp_0[:,t], real.(all_circ_0[t]), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
end 
lines!(ax1, [0, 30], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
scatter!(ax1_2, Temp_rich, [sum(real.(all_circ_0[t]) .<0)./length(path_0) for t in 1:num_temps], color = ("#F8BA17"), label = "", marker = :diamond, markersize = 15, alpha = 1.0)
Label(f[0,1, TopLeft()], "(a)")
linkxaxes!(ax1,ax1_2); hidespines!(ax1_2)
s1 = [MarkerElement(color = ("#285C93", 0.3), markersize = 10, marker = :circle)]
s2 = [MarkerElement(color = ("#F8BA17", 1), markersize = 15, marker = :diamond)]
# l3 = LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)
Legend(f[0, 1], [s1,s2], tellheight = false, tellwidth = false, ["Real λᴶₘ", "p(Stability)"], halign = :right, valign = :center, orientation = :horizontal, framevisible = false, labelsize = 25)
ax2 = Axis(f[1,2], xlabel = "Temperature (°C)", ylabel = "λᴴₘ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2_2 = Axis(f[1,2], ylabel = "p(Reactivity)", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
for t in 1:num_temps
    scatter!(ax2, temp_0[:,t], real.(all_leadH_0[t]), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
end
lines!(ax2, [0, 30], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
scatter!(ax2_2, Temp_rich, [sum(all_leadH_0[t] .>0)./length(path_0) for t in 1:num_temps], color = ("#F8BA17"), label = "", marker = :diamond, markersize = 15, alpha = 1.0)
linkxaxes!(ax2,ax2_2); hidespines!(ax2_2)
Label(f[0,2, TopLeft()], "(b)")
s1 = [MarkerElement(color = ("#285C93", 0.3), markersize = 10, marker = :circle)]
s2 = [MarkerElement(color = ("#F8BA17", 1), markersize = 15, marker = :diamond)]
# l3 = LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)
Legend(f[0, 2], [s1,s2], tellheight = false, tellwidth = false, ["λᴴₘ", "p(Reactivity)"], halign = :right, valign = :center, orientation = :horizontal, framevisible = false, labelsize = 25)
ax3 = Axis(f[1,3], xlabel = "Real λᴶₘ ", ylabel = "λᴴₘ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
scatter!(ax3, real.(vcat(all_circ_0...)), vcat(all_leadH_0...), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
lines!(ax3, [ minimum(real.(vcat(all_circ_0...))), maximum(real.(vcat(all_circ_0...)))], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
lines!(ax3, [0,0], [minimum( vcat(all_leadH_0...)), maximum( vcat(all_leadH_0...))], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
# rea_sta = count(x -> x[1] < 0 && x[2] > 0, zip(real.(vcat(all_circ_0...)), vcat(all_leadH_0...)))/ length(circ_vec)
# text!(ax3, minimum(vcat(all_leadH_0...)), 1.0, text = "$(round(rea_sta * 100, digits = 2))%", align = (:left, :center),fontsize = 20)
Label(f[0,3, TopLeft()], "(c)")
f
save("../results/sta_react_0.pdf", f) 

f = Figure(fontsize = 30, size = (1800, 600));
Label(f[:,0], "Maximal Trade-off", fontsize = 50, rotation = pi/2)
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Real λᴶₘ ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax1_2 = Axis(f[1,1], ylabel = "p(Stability)", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
for t in 1: num_temps
    scatter!(ax1, temp_1[:,t], real.(all_circ_1[t]), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
end 
lines!(ax1, [0, 30], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
scatter!(ax1_2, Temp_rich, [sum(real.(all_circ_1[t]) .<0)./length(path_1) for t in 1:num_temps], color = ("#F8BA17"), label = "", marker = :diamond, markersize = 15, alpha = 1.0)
Label(f[0,1, TopLeft()], "(d)")
linkxaxes!(ax1,ax1_2); hidespines!(ax1_2)
s1 = [MarkerElement(color = ("#285C93", 0.3), markersize = 10, marker = :circle)]
s2 = [MarkerElement(color = ("#F8BA17", 1), markersize = 15, marker = :diamond)]
# l3 = LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)
Legend(f[0, 1], [s1,s2], tellheight = false, tellwidth = false, ["Real λᴶₘ", "p(Stability)"], halign = :right, valign = :center, orientation = :horizontal, framevisible = false, labelsize = 25)
ax2 = Axis(f[1,2], xlabel = "Temperature (°C)", ylabel = "λᴴₘ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2_2 = Axis(f[1,2], ylabel = "p(Reactivity)", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
for t in 1:num_temps
    scatter!(ax2, temp_1[:,t], real.(all_leadH_1[t]), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
end
lines!(ax2, [0, 30], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
scatter!(ax2_2, Temp_rich, [sum(all_leadH_1[t] .>0)./length(path_1) for t in 1:num_temps], color = ("#F8BA17"), label = "", marker = :diamond, markersize = 15, alpha = 1.0)
linkxaxes!(ax2,ax2_2); hidespines!(ax2_2)
Label(f[0,2, TopLeft()], "(e)")
s1 = [MarkerElement(color = ("#285C93", 0.3), markersize = 10, marker = :circle)]
s2 = [MarkerElement(color = ("#F8BA17", 1), markersize = 15, marker = :diamond)]
# l3 = LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)
Legend(f[0, 2], [s1,s2], tellheight = false, tellwidth = false, ["λᴴₘ", "p(Reactivity)"], halign = :right, valign = :center, orientation = :horizontal, framevisible = false, labelsize = 25)
ax3 = Axis(f[1,3], xlabel = "Real λᴶₘ ", ylabel = "λᴴₘ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
scatter!(ax3, real.(vcat(all_circ_1...)), vcat(all_leadH_1...), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
lines!(ax3, [ minimum(real.(vcat(all_circ_1...))), maximum(real.(vcat(all_circ_1...)))], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
lines!(ax3, [0,0], [minimum( vcat(all_leadH_1...)), maximum( vcat(all_leadH_1...))], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
Label(f[0,3, TopLeft()], "(f)")
f
save("../results/sta_react_1.pdf", f) 


######## Community level CUE ###########

path = glob("Eff_iters*", "../data/Eff_p_re_ri/")
all_com_CUE_collect_re = Vector{Vector{Float64}}()
for j in 1: num_temps
    all_com_CUE_H = Float64[];
        for i in 1:length(path)
            @load path[i] all_com_CUE
            push!(all_com_CUE_H, all_com_CUE[j])
        end 
    push!(all_com_CUE_collect_re, all_com_CUE_H)
end 
com_CUE_0 = [mean(all_com_CUE_collect_0[t]) for t in 1:num_temps]
com_CUE_err_0 = [std(all_com_CUE_collect_0[t])/sqrt(length(all_com_CUE_collect_0[t])) for t in 1: num_temps]
com_CUE_1 = [mean(all_com_CUE_collect_1[t]) for t in 1:num_temps]
com_CUE_err_1 = [std(all_com_CUE_collect_1[t])/sqrt(length(all_com_CUE_collect_1[t])) for t in 1: num_temps]
com_CUE_re = [mean(all_com_CUE_collect_re[t]) for t in 1:num_temps]
com_CUE_err_re = [std(all_com_CUE_collect_re[t])/sqrt(length(all_com_CUE_collect_re[t])) for t in 1: num_temps]


f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Community-level CUE", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, com_CUE_0, color = ("#EF8F8C",0.8), linewidth = 5, label = "ρ = 0")
band!(ax1, Temp_rich, com_CUE_0 .- com_CUE_err_0, com_CUE_0 .+ com_CUE_err_0, color = ("#EF8F8C", 0.3))
lines!(ax1, Temp_rich, com_CUE_re, color = ("#376298",0.8), linewidth = 5, label = "Realistic ρ")
band!(ax1, Temp_rich, com_CUE_re .- com_CUE_err_re, com_CUE_0 .+ com_CUE_err_re, color = ("#376298", 0.3))
lines!(ax1, Temp_rich, com_CUE_1, color = ("#4F363E",0.8), linewidth = 5, label = "ρ = -1")
band!(ax1, Temp_rich, com_CUE_1 .- com_CUE_err_1, com_CUE_1 .+ com_CUE_err_1, color = ("#4F363E", 0.3))
axislegend(position = :rt)
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/com_CUE.pdf", f) 

###################################
path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_ii_collect = Vector{Vector{Float64}}(); all_ij_collect = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_ii_H = Float64[]; all_ij_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij #all_ℵii_sur
        append!(all_ii_H, all_ℵii[j]); append!(all_ij_H, all_ℵij[j])
        next!(progress)
    end 
    push!(all_ii_collect, all_ii_H); push!(all_ij_collect, all_ij_H)
end 

temp_ii = [repeat([Temp_rich[t]], length(all_ii_collect[t])) for t in 1:31]
temp_ij = [repeat([Temp_rich[t]], length(all_ij_collect[t])) for t in 1:31]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature ", ylabel = "α", xlabelsize = 50, ylabelsize = 50)
# for i in 1: 31
#     scatter!(ax, real.(all_circ[i]), real.(all_leadH[i]), color = cs[i], markersize = 10, alpha = 0.05)
# end 
boxplot!(ax, vcat(temp_ii...), log.(abs.(vcat(all_ii_collect...))), color = ("#FA8328", 0.7), label = "αᵢᵢ") 
boxplot!(ax, vcat(temp_ij...), log.(abs.(vcat(all_ij_collect...))), color = ("#069F66", 0.7), label = "αᵢⱼ") 
# lines!(ax, [ minimum(real.(vcat(all_circ...))), maximum(real.(vcat(all_circ...)))], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
# lines!(ax, [0,0], [minimum( vcat(all_leadH...)), maximum( vcat(all_leadH...))], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
# # Colorbar(f[1,2], colorrange = [0, num_temps], colormap = cs, label = "Temperature")
# Label(f[1,1, TopLeft()], "(b)")
axislegend(position = :lt)
f
save("../results/ijii0.pdf", f) 


RO0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 29]
RO0_err = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 30]
RO1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 29]
RO1_err = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 29]


f = Figure(fontsize = 35, size = (1500, 900));
ax1 = Axis(f[1,1], title = "Minimum Overlap", xlabel = "α", ylabel = "frequency (αii)", xlabelsize = 25, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "", xlabelsize = 25, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
hist!(ax1, all_ii_collect[argmin(RO0)], color = ("#FA8328", 0.7), bins = 100)
hist!(ax2, all_ij_collect[argmin(RO0)], color = ("#069F66", 0.7), bins = 100)
linkxaxes!(ax1,ax2)
ax3 = Axis(f[1,2], title = "Maximum Overlap", xlabel = "α", ylabel = "", xlabelsize = 25, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax4 = Axis(f[1,2], ylabel = "", xlabelsize = 25, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax4)
hist!(ax3, all_ii_collect[argmax(RO0)], color = ("#FA8328", 0.7), bins = 100)
hist!(ax4, all_ij_collect[argmax(RO0)], color = ("#069F66", 0.7), bins = 100)
linkxaxes!(ax3,ax4)
p1 = [PolyElement(color = ("#FA8328", 0.7), strokecolor = :transparent)]
p2 = [PolyElement(color = ("#069F66", 0.7), strokecolor = :transparent)]
Legend(f[1,1], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
Legend(f[1,2], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
f

std(all_ii_collect[argmin(RO0)])
std(all_ij_collect[argmin(RO0)])
std(all_ii_collect[argmax(RO0)])
std(all_ij_collect[argmax(RO0)])

#########################
path = glob("Eff_iters*", "../data/Eff_p0_cos/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_RO_collect = Vector{Vector{Float64}}(); all_ulO_collect = Vector{Vector{Float64}}(); all_Rul_collect = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_RO_H = Float64[]; all_ulO_H = Float64[]; all_Rul_H = Float64[]
    for i in 1:length(path)
        @load path[i] RO ulO Rul 
        append!(all_RO_H, RO[j]); append!(all_ulO_H, ulO[j]); append!(all_Rul_H, Rul[j])
        next!(progress)
    end 
    push!(all_RO_collect, all_RO_H); push!(all_ulO_collect, all_ulO_H); push!(all_Rul_collect, all_Rul_H)
end 
RO_0 = [mean(all_RO_collect[t]) for t in 1:num_temps]
RO_err_0 = [std(all_RO_collect[t])/sqrt(length(all_RO_collect[t])) for t in 1: num_temps]
ulO_0 = [mean(all_ulO_collect[t]) for t in 1:num_temps]

path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_Eu_collect = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_Eu_H = Float64[]
    for i in 1:length(path)
        @load path[i] all_Eu
        pair_E = [all_Eu[j][i] + all_Eu[j][s] for i in 1:N for s in 1:N]
        append!(all_Eu_H, pair_E)
        next!(progress)
    end 
    push!(all_Eu_collect, all_Eu_H)
end 

hist(vcat(all_Eu_collect...))
all_Eu_c = vcat(all_Eu_collect...)

########### 
path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
dℵij_collect = Vector{Vector{Float64}}() ; all_R_collect = Vector{Vector{Float64}}()
# all_αijii_collect = Vector{Vector{Float64}}()
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    dℵij_H = Float64[]; all_R_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij all_R
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(dℵij_H, dℵij); append!(all_R_H, all_R[j])
        next!(progress)
    end 
    push!(dℵij_collect, dℵij_H); push!(all_R_collect, all_R_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

αij_d = [mean(dℵij_collect[t]) for t in 1:num_temps]
αij_d_err = [std(dℵij_collect[t])/sqrt(length(dℵij_collect[t])) for t in 1: num_temps]
meanR = [mean(all_R_collect[t]) for t in 1:num_temps]
R_err = [std(all_R_collect[t])/sqrt(length(all_R_collect[t])) for t in 1: num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
# Label(f[0,1], "Minimal Trade-off", fontsize = 50)
ax1 = Axis(f[1,1], title = "Minimal Trade-off", xlabel = "Temperature (°C)", ylabel = "Effective Resource Competition\n(αji/αii)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Resource Abundance", yaxisposition = :right, yticklabelalign = (:left, :center), xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
# hidedecorations!(ax3, grid = false, ticks = true, ticklabels = true)
lines!(ax1, Temp_rich, αij_d, color = ("#376298",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, αij_d .- αij_d_err , αij_d.+ αij_d_err, color = ("#376298", 0.3))
lines!(ax2, Temp_rich, meanR, color =( "#F8BA17", 0.9), linewidth = 5, label = "")
band!(ax2, Temp_rich, meanR .- R_err , meanR .+ R_err , color = ("#F8BA17", 0.5))
linkxaxes!(ax1,ax2)
lines!(ax1, [0, 30], [1, 1], linestyle = :dash, color = ("#4F363E", 0.9), linewidth = 2)
text!(ax1, 0, 1.05, text = "↑ αⱼᵢ > αᵢᵢ", align = (:left, :center),fontsize = 30)
text!(ax1, 0, 0.95, text = "↓ αⱼᵢ < αᵢᵢ", align = (:left, :center),fontsize = 30)
l1 = [LineElement(color = ("#376298", 0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#F8BA17", 0.9), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αji/αii", "Resource"], halign = :center, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/αR0.pdf", f) 


############### ρ ################
path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
# all_ρ_collect = Vector{Vector{Float64}}(); all_V_collect = Vector{Vector{Float64}}(); all_E_collect = Vector{Vector{Float64}}();
# all_upper_collect = Vector{Vector{Float64}}(); all_lower_collect = Vector{Vector{Float64}}(); 
all_var_d_collect_1 = Vector{Vector{Float64}}(); all_var_offd_collect_1 = Vector{Vector{Float64}}(); 
all_May_1 = Vector{Vector{Float64}}(); all_Tang_Allesina_1 = Vector{Vector{Float64}}();
all_A_diag_1 = Vector{Vector{Float64}}(); all_A_off_1 = Vector{Vector{Float64}}();
# all_A_diag_sur_1 = Vector{Vector{Float64}}(); all_A_off_sur_1 = Vector{Vector{Float64}}();
all_Jac_diag_1 = Vector{Vector{Float64}}(); all_Jac_off_1 = Vector{Vector{Float64}}()

@time for j in 1: num_temps
    # all_ρ_H = Float64[]; all_V_H = Float64[]; all_E_H = Float64[]; all_upper_H =  Float64[]; all_lower_H = Float64[]; 
    all_var_d_H = Float64[]; all_var_offd_H = Float64[]; all_May_H = Float64[]; all_Tang_Allesina_H = Float64[];
    A_diag_H = Float64[]; A_off_H = Float64[]; # A_diag_sur_H = Float64[]; A_off_sur_H = Float64[]; 
    Jac_diag_H = Float64[];Jac_off_H = Float64[]

    for i in 1:length(path)
        @load path[i] all_ℵii all_ℵij all_r all_ℵii_sur all_C
        sur = findall(x -> x in all_ℵii_sur[j], all_ℵii[j])
        N_s = length(sur)
        C = all_C[j]; r = all_r[j]
        C_all = zeros(N)
        C_all[sur] = C
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        A_sur = A[sur, sur]
        # A_off_sur = [A_sur[i, j] for i in 1:N_s for j in 1:N_s if j != i]
        # LV_Jac = [A[i, j]*C_all[i] for i in 1:N, j in 1:N]
        # LV_Jac[diagind(LV_Jac)] .= [r[i] + A[i, i]*C_all[i] + sum(A[i, j]*C_all[j] for j in 1:N) for i in 1:N]
        LV_Jac = [A_sur[i, j]*C[i] for i in 1:N_s, j in 1:N_s]
        off_diag_Jac = [LV_Jac[i, j] for i in 1:N_s for j in 1:N_s if j != i]
        upper = [LV_Jac[i, j] for i in 1:N_s for j in 1:N_s if j > i]
        lower = [LV_Jac[j, i] for i in 1:N_s for j in 1:N_s if j > i]
        E = mean(off_diag_Jac); V = var(off_diag_Jac)
        # E2 = mean(upper .* lower)
        # ρ = (E2-E^2)/V
        d = mean(diag(LV_Jac))
        θ = d/V
        May = std(off_diag_Jac)*sqrt(N_s) - d # < 0
        Tang_Allesina = sqrt(N_s) - θ/ (1 + E^2/V^2) # < 0
        # append!(all_ρ_H,  ρ); append!(all_V_H, V);append!(all_E_H, E); append!(all_upper_H, upper); append!(all_lower_H, lower); 
        append!(all_var_d_H, std(diag(LV_Jac))); append!(all_var_offd_H, std(off_diag_Jac))
        append!(all_May_H, May); append!(all_Tang_Allesina_H, Tang_Allesina);
        append!(A_diag_H, std(all_ℵii[j])); append!(A_off_H, std(all_ℵij[j])); 
        # append!(A_diag_sur_H, diag(A_sur)); append!(A_off_sur_H, A_off_sur); 
        append!(Jac_diag_H, std(diag(LV_Jac))); append!(Jac_off_H, std(off_diag_Jac))
        next!(progress)
    end 
    # push!(all_ρ_collect, all_ρ_H); push!(all_V_collect, all_V_H);push!(all_E_collect, all_E_H); push!(all_upper_collect, all_upper_H); push!(all_lower_collect, all_lower_H);
    push!(all_var_d_collect_1, all_var_d_H); push!(all_var_offd_collect_1, all_var_offd_H);
    push!(all_May_1, all_May_H); push!(all_Tang_Allesina_1, all_Tang_Allesina_H);
    # push!(all_A_diag_sur_1, A_diag_sur_H); push!(all_A_off_sur_1, A_off_sur_H);
    push!(all_A_diag_1, A_diag_H); push!(all_A_off_1, A_off_H);
    push!(all_Jac_diag_1, Jac_diag_H); push!(all_Jac_off_1, Jac_off_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

var_Jij_0 = [mean(all_var_offd_collect_0[t]) for t in 1: num_temps]
var_Jij_err_0 = [std(all_var_offd_collect_0[t])/sqrt(length(all_var_offd_collect_0[t])) for t in 1: num_temps]
var_Jii_0 = [mean(all_var_d_collect_0[t]) for t in 1: num_temps]
var_Jii_err_0 = [std(all_var_d_collect_0[t])/sqrt(length(all_var_d_collect_0[t])) for t in 1: num_temps]

var_Jijii_0 = [mean(all_var_offd_collect_0[t] - all_var_d_collect_0[t]) for t in 1: num_temps]
var_Jijii_err_0 = [std(all_var_offd_collect_0[t] - all_var_d_collect_0[t])/sqrt(length(all_var_offd_collect_0[t])) for t in 1: num_temps]
var_Jijii_1 = [mean(all_var_offd_collect_1[t] - all_var_d_collect_1[t]) for t in 1: num_temps]
var_Jijii_err_1 = [std(all_var_offd_collect_1[t] - all_var_d_collect_1[t])/sqrt(length(all_var_offd_collect_1[t])) for t in 1: num_temps]

var_Jij_1 = [mean(all_var_offd_collect_1[t]) for t in 1: num_temps]
var_Jij_err_1 = [std(all_var_offd_collect_1[t])/sqrt(length(all_var_offd_collect_1[t])) for t in 1: num_temps]
var_Jii_1 = [mean(all_var_d_collect_1[t]) for t in 1: num_temps]
var_Jii_err_1 = [std(all_var_d_collect_1[t])/sqrt(length(all_var_d_collect_1[t])) for t in 1: num_temps]

var_Jii = vcat(var_Jii_0, var_Jii_1)
var_Jij = vcat(var_Jij_0, var_Jij_1)
var_Jijii = vcat(var_Jijii_0, var_Jijii_1)
var_Jii_err = vcat(var_Jii_err_0, var_Jii_err_1)
var_Jij_err = vcat(var_Jij_err_0, var_Jij_err_1)
var_Jijii_err = vcat(var_Jijii_err_0, var_Jijii_err_1)

May_0 = [mean(all_May_0[t]) for t in 1: num_temps]
May_err_0 = [std(all_May_0[t])/sqrt(length(all_May_0[t])) for t in 1: num_temps]
May_1 = [mean(all_May_1[t]) for t in 1: num_temps]
May_err_1 = [std(all_May_1[t])/sqrt(length(all_May_1[t])) for t in 1: num_temps]
May = vcat(May_0, May_1)
May_err = vcat(May_err_0, May_err_1)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1],title = "May's Criterion (stable if <0)", xlabel = L"σ_{J}\sqrt{Sc} - \bar{d}", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, May_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(May_0, stability_0, May_err_0)
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    cap_length_0 = 0.001 * mean(stability_0)  # Length of horizontal caps
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, May_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(May_1, stability_1, May_err_1)
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    cap_length_0 = 0.001 * mean(stability_1)  # Length of horizontal caps
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rb)
f
save("../results/May_Criteria.pdf", f) 

Tang_Allesina_0 = [mean(all_Tang_Allesina_0[t]) for t in 1: num_temps]
Tang_Allesina_err_0 = [std(all_Tang_Allesina_0[t])/sqrt(length(all_Tang_Allesina_0[t])) for t in 1: num_temps]
Tang_Allesina_1 = [mean(all_Tang_Allesina_1[t]) for t in 1: num_temps]
Tang_Allesina_err_1 = [std(all_Tang_Allesina_1[t])/sqrt(length(all_Tang_Allesina_1[t])) for t in 1: num_temps]
Tang_Allesina = vcat(Tang_Allesina_0, Tang_Allesina_1)
Tang_Allesina_err = vcat(Tang_Allesina_err_0, Tang_Allesina_err_1)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1],title = "Allesina & Tang's Criterion (stable if <0)", xlabel = L"\sqrt{Sc} - θ/(1+E^2(J)/σ_{J}^2)", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, Tang_Allesina_0, stability_0, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(Tang_Allesina_0, stability_0, Tang_Allesina_err_0)
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    cap_length_0 = 0.001 * mean(stability_0)  # Length of horizontal caps
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
scatter!(ax, Tang_Allesina_1, stability_1, color = "#9A2B1A", markersize = 15, alpha = 0.8, label = "ρ = -1")
for (x, y, e) in zip(Tang_Allesina_1, stability_1, Tang_Allesina_err_1)
    lines!(ax, [x - e, x + e], [y, y], color = ("#9A2B1A", 0.4), linewidth = 1)
    cap_length_0 = 0.001 * mean(stability_1)  # Length of horizontal caps
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#9A2B1A", 0.4), linewidth = 1)
end
axislegend(position = :rb)
f
save("../results/Tang_Allesina_Criteria.pdf", f) 

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = L"σ(J_{ij}) - σ(J_{ii})", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, var_Jijii, stability_all, color = "#4F363E", markersize = 15, alpha = 0.8)
for (x, y, e) in zip(var_Jijii, stability_all, var_Jijii_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#4F363E", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#4F363E", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#4F363E", 0.4), linewidth = 1)
end
# axislegend(position = :rc)
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/var_Jijii_sta.pdf", f) 


std_diag_J_0 = [mean(all_Jac_diag_0[t]) for t in 1: num_temps]
std_diag_J_err_0 = [std(all_Jac_diag_0[t])/sqrt(length(all_Jac_diag_0[t])) for t in 1: num_temps]
std_off_J_0 = [mean(all_Jac_off_0[t]) for t in 1: num_temps]
std_off_J_err_0 = [std(all_Jac_off_0[t])/sqrt(length(all_Jac_off_0[t])) for t in 1: num_temps]
std_diag_A_0 = [mean(all_A_diag_0[t]) for t in 1: num_temps]
std_diag_A_err_0 = [std(all_A_diag_0[t])/sqrt(length(all_A_diag_0[t])) for t in 1: num_temps]
std_off_A_0 = [mean(all_A_off_0[t]) for t in 1: num_temps]
std_off_A_err_0 = [std(all_A_off_0[t])/sqrt(length(all_A_off_0[t])) for t in 1: num_temps]

std_diag_J_1 = [mean(all_Jac_diag_1[t]) for t in 1: num_temps]
std_diag_J_err_1 = [std(all_Jac_diag_1[t])/sqrt(length(all_Jac_diag_1[t])) for t in 1: num_temps]
std_off_J_1 = [mean(all_Jac_off_1[t]) for t in 1: num_temps]
std_off_J_err_1 = [std(all_Jac_off_1[t])/sqrt(length(all_Jac_off_1[t])) for t in 1: num_temps]
std_diag_A_1 = [mean(all_A_diag_1[t]) for t in 1: num_temps]
std_diag_A_err_1 = [std(all_A_diag_1[t])/sqrt(length(all_A_diag_1[t])) for t in 1: num_temps]
std_off_A_1 = [mean(all_A_off_1[t]) for t in 1: num_temps]
std_off_A_err_1 = [std(all_A_off_1[t])/sqrt(length(all_A_off_1[t])) for t in 1: num_temps]

std_diag_J = vcat(std_diag_J_0, std_diag_J_1)
std_diag_J_err = vcat(std_diag_J_err_0, std_diag_J_err_1)
std_off_J = vcat(std_off_J_0, std_off_J_1)
std_off_J_err = vcat(std_off_J_err_0, std_off_J_err_1)
std_diag_A = vcat(std_diag_A_0, std_diag_A_1)
std_diag_A_err = vcat(std_diag_A_err_0, std_diag_A_err_1)
std_off_A = vcat(std_off_A_0, std_off_A_1)
std_off_A_err = vcat(std_off_A_err_0, std_off_A_err_1)

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = L"σ(A_{ii})", ylabel = L"σ(J_{feas})", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, std_diag_A, std_diag_J, color = "#FA8328", markersize = 15, alpha = 0.8, label = "diagonal")
for (x, y, e_A, e_J) in zip(std_diag_A, std_diag_J, std_diag_A_err, std_diag_J_err)
    # Vertical line
    lines!(ax, [x - e_A, x + e_A], [y, y], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x, x], [y - e_J, y + e_J], color = ("#FA8328", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_y = 0.02 * mean(std_diag_J)  # Length of horizontal caps
    cap_length_x = 0.02 * mean(std_diag_A)
    lines!(ax, [x - e_A, x - e_A], [y - cap_length_y, y + cap_length_y], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x + e_A, x + e_A], [y - cap_length_y, y + cap_length_y], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_x, x + cap_length_x], [y - e_J, y - e_J], color = ("#FA8328", 0.4), linewidth = 1)
    lines!(ax, [x - cap_length_x, x + cap_length_x], [y + e_J, y + e_J], color = ("#FA8328", 0.4), linewidth = 1)
end
ax1 = Axis(f[1,1], xlabel = L"σ(A_{ij})", xlabelsize = 50, ylabelsize = 50, xaxisposition = :top, ygridvisible = false, xgridvisible = false,yticklabelsvisible = false, ylabelvisible = false)
scatter!(ax1, std_off_A, std_off_J, color = "#069F66", markersize = 15, alpha = 0.8, label = "off-diagonal")
for (x, y, e_A, e_J) in zip(std_off_A, std_off_J, std_off_A_err, std_off_J_err)
    # Vertical line
    lines!(ax1, [x - e_A, x + e_A], [y, y], color = ("#069F66", 0.4), linewidth = 1)
    lines!(ax1, [x, x], [y - e_J, y + e_J], color = ("#069F66", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_y = 0.02 * mean(std_off_J)  # Length of horizontal caps
    cap_length_x = 0.02 * mean(std_off_A)
    lines!(ax1, [x - e_A, x - e_A], [y - cap_length_y, y + cap_length_y], color = ("#069F66", 0.4), linewidth = 1)
    lines!(ax1, [x + e_A, x + e_A], [y - cap_length_y, y + cap_length_y], color = ("#069F66", 0.4), linewidth = 1)
    lines!(ax1, [x - cap_length_x, x + cap_length_x], [y - e_J, y - e_J], color = ("#069F66", 0.4), linewidth = 1)
    lines!(ax1, [x - cap_length_x, x + cap_length_x], [y + e_J, y + e_J], color = ("#069F66", 0.4), linewidth = 1)
end
linkyaxes!(ax,ax1); hidespines!(ax1)
s1 = [MarkerElement(color = ("#FA8328", 0.8), markersize = 15, marker = :circle)]
s2 = [MarkerElement(color = ("#069F66", 0.8), markersize = 15, marker = :circle)]
Legend(f[1,1], [s1,s2], tellheight = false, tellwidth = false, ["diagonal", "off-diagonal"], halign = :right, valign = :bottom)
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/σJ_feas_σA.pdf", f) 


plot([mean(all_V_collect[t] .* 10 ) for t in 1:num_temps])

up = [mean(all_upper_collect[t]) for t in 1:num_temps]
up_err = [std(all_upper_collect[t])/sqrt(length(all_upper_collect[t])) for t in 1:num_temps]
low = [mean(all_lower_collect[t]) for t in 1:num_temps]
low_err = [std(all_lower_collect[t])/sqrt(length(all_lower_collect[t])) for t in 1:num_temps]
plot(all_upper_collect[1], all_lower_collect[1])

scatter(vcat(all_ρ_collect...))
ρ_ijji_t = [mean(all_ρ_collect[t]) for t in 1:num_temps]
ρ_ijji_err = [std(all_ρ_collect[t])/sqrt(length(all_ρ_collect[t])) for t in 1:num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "ρ(α)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, ρ_ijji_t, color = ("#376298",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, ρ_ijji_t .- ρ_ijji_err , ρ_ijji_t.+ ρ_ijji_err, color = ("#376298", 0.3))
# Label(f[1,1, TopLeft()], "(b)")
f

stability = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 59]
# stability = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 59]

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "ρ(α)", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, ρ_ijji_t, stability, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
for (x, y, e) in zip(ρ_ijji_t, stability, ρ_ijji_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#376298", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color = ("#376298", 0.4), linewidth = 1)
end
axislegend(position = :lb)
# Label(f[1,1, TopLeft()], "(a)")
f



all_ii_collect = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_ii = Float64[]
    for i in 1:length(path)
        @load path[i] all_ℵii
        mean_ℵii = mean(all_ℵii[j])
        append!(all_ii, mean_ℵii)
    end 
    push!(all_ii_collect, all_ii)
end

# diff_t = [mean(sqrt.(N .* all_V_collect[t]) .* (1 .+ all_ρ_collect[t]) .- all_E_collect[t]) for t in 1:num_temps]
# diff_d_t = [mean(sqrt.(N .* all_V_collect[t]) .* (1 .+ all_ρ_collect[t]) .- all_E_collect[t] .- (-all_ii_collect[t])) for t in 1:num_temps]

# vcat(all_ρ_collect...)
# diff_d = sqrt.(N .* vcat(all_V_collect...)) .* (1 .+ vcat(all_ρ_collect...)) .- vcat(all_E_collect...)
# vcat(all_E_collect...)

# f = Figure(fontsize = 35, size = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Jᵢⱼ", ylabel = "Jⱼᵢ", xlabelsize = 50, ylabelsize = 50)
# scatter!(ax, up, low, color = "#376298", markersize = 15, alpha = 0.8, label = "ρ = 0")
# for (x, y, e_u, e_l) in zip(up, low, up_err, low_err)
#     # Vertical line
#     lines!(ax, [x - e_u, x + e_u], [y, y], color = ("#376298", 0.4), linewidth = 1)
#     lines!(ax, [x, x], [y - e_l, y + e_l], color = ("#376298", 0.4), linewidth = 1)
#     # Horizontal caps
#     cap_length = 0.001 * mean(up)  # Length of horizontal caps
#     # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
#     # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
#     lines!(ax, [x - e_u, x - e_u], [y - cap_length, y + cap_length], color = ("#376298", 0.4), linewidth = 1)
#     lines!(ax, [x + e_u, x + e_u], [y - cap_length, y + cap_length], color = ("#376298", 0.4), linewidth = 1)
#     lines!(ax, [x - cap_length, x + cap_length], [y - e_l, y - e_l], color = ("#376298", 0.4), linewidth = 1)
#     lines!(ax, [x - cap_length, x + cap_length], [y + e_l, y + e_l], color = ("#376298", 0.4), linewidth = 1)
# end
# # axislegend(position = :rb)
# # Label(f[1,1, TopLeft()], "(a)")
# f

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


