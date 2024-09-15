include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob
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
CSV.write("../results/Eff_results_p-1_new.csv", Eff_results, writeheader=false)
@save "../results/Feas_CR_dist_p-1_new.jld2" all_Rrela_collect all_Crela_collect all_R_collect all_C_collect all_ii_collect all_ij_collect all_ii_sur_collect all_ij_sur_collect all_r_collect all_r_sur_collect

##### Loading #######
# Eff_results = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)
# @load "../results/Feas_CR_dist_p0_new.jld2" all_Rrela_collect all_Crela_collect all_R_collect all_C_collect all_ii_collect all_ij_collect all_ii_sur_collect all_ij_sur_collect all_r_collect all_r_sur_collect

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
path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_rich_collect_1 = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    all_rich_H = Float64[];
    # all_αijii_H = Float64[]
    for i in 1:length(path)
        @load path[i] all_ℵii_sur
        append!(all_rich_H, length(all_ℵii_sur[j]))
        # append!(all_αijii_H,  αijii);
        next!(progress)
    end 
    push!(all_rich_collect_1, all_rich_H)
    # push!(all_αijii_collect, all_αijii_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

stability_0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 59]
stability_1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 59]
mean_rich_0 = [mean(all_rich_collect_0[t]) for t in 1:num_temps]
rich_err_0 = [std(all_rich_collect_0[t])/sqrt(length(all_rich_collect_0[t])) for t in 1: num_temps]
mean_rich_1 = [mean(all_rich_collect_1[t]) for t in 1:num_temps]
rich_err_1 = [std(all_rich_collect_1[t])/sqrt(length(all_rich_collect_1[t])) for t in 1: num_temps]

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
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/rich_sta.pdf", f) 

################# evenness ##################
function pielou_evenness(C_rela)
    H = -sum(C_rela .* log.(C_rela))
    H_max = log(length(C_rela))
    return H / H_max
end

path = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_even_collect_1 = Vector{Vector{Float64}}()
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
    push!(all_even_collect_1, all_even_H)
    # push!(all_αijii_collect, all_αijii_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

even_mean = [mean(all_even_collect[t]) for t in 1: num_temps]
even_err = [std(all_even_collect[t])/sqrt(length(all_even_collect[t])) for t in 1: num_temps]
Eff_results.stability
f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Evenness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, even_mean, Eff_results.stability, color = "#285C93", markersize = 15, alpha = 0.8)
for (x, y, e) in zip(even_mean, Eff_results.stability, even_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = ("#285C93", 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(Eff_results.stability)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = ("#285C93", 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = ("#285C93", 0.4), linewidth = 1)
end
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/sta_varC-1.pdf", f) 


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
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/sta_varC.pdf", f) 

############### σ vs. stability 
stability_0 = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 59]
stability_1 = CSV.read("../results/Eff_results_p-1_new.csv", DataFrame, header=false)[:, 59]

path = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
all_ij_collect_0 = Vector{Vector{Float64}}()
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_ij_H = Float64[]
    for i in 1:length(path)
        @load path[i]  all_ℵii all_ℵij #all_ℵii_sur
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        ij = std(log.(abs.([A[i, j]/A[i, i] for i in 1:N for j in 1:N if j != i])))
        append!(all_ij_H, ij)
        next!(progress)
    end 
    push!(all_ij_collect_0, all_ij_H)
end 

var_ij_0 = [mean(all_ij_collect_0[t]) for t in 1: num_temps]
var_err_0 = [std(all_ij_collect_0[t])/sqrt(length(all_ij_collect_0[t])) for t in 1: num_temps]
var_ij_1 = [mean(all_ij_collect_1[t]) for t in 1: num_temps]
var_err_1 = [std(all_ij_collect_1[t])/sqrt(length(all_ij_collect_1[t])) for t in 1: num_temps]

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
save("../results/sta_σ.pdf", f) 


########### upper lower #############

path = glob("Eff_iters*", "../data/Eff_p0_new/")
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
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/sta_diag_dom.pdf", f) 

##########

path_0 = glob("Eff_iters*", "../data/Eff_p0_new/")
progress = Progress(length(path_0)*num_temps; desc="Progress running:")
num_temps = 31
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
all_leading_collect0 = Vector{Vector{ComplexF64}}()
@time for j in 1: num_temps
    all_leading_H = ComplexF64[]
    for i in 1:length(path_0)
        @load path_0[i] all_ℵii all_ℵij all_r_sur all_ℵii_sur all_C
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
        append!(all_leading_H,  leading)
        next!(progress)
    end 
    push!(all_leading_collect0, all_leading_H)
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
path = glob("Eff_iters*", "../data/Eff_p-1_new/")


all_circ = Vector{Vector{ComplexF64}}()
for j in 1: num_temps
    # if (j-1) % 5 == 0
        circ_leading_H = ComplexF64[]
        for i in 1:length(path)
            @load path[i] all_leading diag_dominance
            push!(circ_leading_H, all_leading[j])
        end 
    # end 
    push!(all_circ, circ_leading_H)
end 

temp = hcat([repeat([Temp_rich[t]], length(all_circ[t])) for t in 1:num_temps]...)


f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Leading Eigen (real)", xlabelsize = 50, ylabelsize = 50)
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



######## Community level CUE ###########

path = glob("Eff_iters*", "../data/Eff_p0_ri/")
all_com_CUE_collect = Vector{Vector{Float64}}()
for j in 1: num_temps
    all_com_CUE_H = Float64[];
        for i in 1:length(path)
            @load path[i] all_com_CUE
            push!(all_com_CUE_H, all_com_CUE[j])
        end 
    push!(all_com_CUE_collect, all_com_CUE_H)
end 

com_CUE_0 = [mean(all_com_CUE_collect[t]) for t in 1:num_temps]
com_CUE_err_0 = [std(all_com_CUE_collect[t])/sqrt(length(all_com_CUE_collect[t])) for t in 1: num_temps]
f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Community-level CUE", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, com_CUE_0, color = ("#EF8F8C",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, com_CUE_0 .- com_CUE_err_0, com_CUE_0 .+ com_CUE_err_0, color = ("#EF8F8C", 0.3))
axislegend(position = :rb)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/com_CUE_0.pdf", f) 

com_CUE_1 = [mean(all_com_CUE_collect[t]) for t in 1:num_temps]
com_CUE_err_1 = [std(all_com_CUE_collect[t])/sqrt(length(all_com_CUE_collect[t])) for t in 1: num_temps]
f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Community-level CUE", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, com_CUE_1, color = ("#EF8F8C",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, com_CUE_1 .- com_CUE_err_1, com_CUE_1 .+ com_CUE_err_1, color = ("#EF8F8C", 0.3))
axislegend(position = :rb)
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/com_CUE-1.pdf", f) 



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
