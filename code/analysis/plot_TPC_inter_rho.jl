include("./fitting.jl");
using ProgressMeter, RCall

##########################################
N=100
M=50
L = fill(0.0, N)
### Temp params 
num_temps = 38

condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

Tr=273.15+10; Ed=3.5 
tspan = (0.0, 2.5e10)
x0 = vcat(fill(0.1, N), fill(1, M)) 
Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)
temp_SS(T, params) = params[1] .* exp.((-params[2]./k) * ((1 ./T) .-(1/Tr)))./(1 .+ (params[2]./(params[4] .- params[2])) .* exp.(params[4]/k * (1 ./params[3] .- 1 ./T)))

tempe_d = reshape(Int.(range(1, N*N, N*N)), N, N)'
tempe_ii = diag(tempe_d)
tempe_ij =[tempe_d[i,j] for i in 1:N for j in 1:N if j!=i]
function filter_fitted(path, p)
    fitted = CSV.read(path, DataFrame, header=false)
    df_names = ["B0","E","Th","Ed","AIC","r2"]
    fitted = rename(fitted, df_names)
    
    if nrow(fitted) == 100  # Diagonal data
        fitted[!, :id] = tempe_ii
        fitted[!, :Bu] = p.B[:,1]
        fitted[!, :Eu] = p.E[:,1]
        fitted[!, :Thu] = p.Tp[:,1]
    elseif nrow(fitted) == 9900  # Off-diagonal data
        fitted[!, :id] = tempe_ij
    else
        @warn "Unexpected number of rows: $(nrow(fitted)). Expected 100 (αii) or 9900 (αij)"
        if nrow(fitted) <= 100
            fitted[!, :id] = tempe_ii[1:nrow(fitted)]
            fitted[!, :Bu] = p.B[1:nrow(fitted),1]
            fitted[!, :Eu] = p.E[1:nrow(fitted),1]
            fitted[!, :Thu] = p.Tp[1:nrow(fitted),1]
        else
            fitted[!, :id] = tempe_ij[1:nrow(fitted)]
        end
    end
    
    # Apply filters
    fitted = fitted[fitted.E .> eps(), :]
    fitted = fitted[fitted.B0 .> eps(), :]
    fitted = fitted[fitted.Th .< 323.15, :]
    fitted = fitted[fitted.Th .> 273.15, :]
    fitted = fitted[fitted.r2 .> 0.9, :]
    
    return fitted
end

idx = collect(CartesianIndices(zeros(Float64, 5, 5)))
ind_diag = diag(idx)
ind_off = [idx[i,j] for i in 1:5 for j in 1:5 if i != j]
all_ind = vcat(ind_diag, ind_off)

p_0 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=273.15 + 25, ρ_t=[0.0000 0.0000], Tr=Tr, Ed=Ed)

##########################################
# progress = Progress(990; desc="Progress running:")
# fitted_all = zeros(Float64, 9900, 6)
# for index in 1: 990
#     next!(progress)
#     # @load "../results/fit_allij_0_new/fit_allij_0_$(index).jld2" fitted
#     @load "../results/inter_ij0_rho/fit_allij_0_rho4_$(index).jld2" fitted
#     i = 10*(index-1) + 10
#     fitted_all[10*(index-1)+1:i,:] = fitted
# end 

# df_names = ["B0","E","Th","Ed","AIC","r2"]
# fitted_all = DataFrame(fitted_all, df_names);
# CSV.write("../results/αij_fitted_rho4.csv", fitted_all, writeheader=false)

# # fitted_all_filter = fitted_all[fitted_all.B0 .> 0, :]
# # fitted_all_filter = fitted_all_filter[fitted_all_filter.E .> 0, :]

####################### Data for plots ##############
### rho2
fitted_ii_0_rho2 = CSV.read("../results/αii_fitted_rho2.csv", DataFrame, header=false)
fitted_ij_0_rho2 = CSV.read("../results/αij_fitted_rho2.csv", DataFrame, header=false)
# fitted = vcat(fitted_ii, fitted_ij)
df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_ii_0_rho2 = DataFrame(fitted_ii_0_rho2, df_names); 
fitted_ij_0_rho2 = DataFrame(fitted_ij_0_rho2, df_names); 

@load "../data/1com0_rho2.jld2" all_ℵii all_ℵij # all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
all_ℵii_0_rho2 = all_ℵii; all_ℵij_0_rho2 = all_ℵij
path_ii_rho2 = ["../results/αii_fitted_rho2.csv"]
path_ij_rho2 = ["../results/αij_fitted_rho2.csv"]
fitted_0ii_rho2 = filter_fitted(path_ii_rho2[1], p_0)
fitted_0ij_rho2 = filter_fitted(path_ij_rho2[1], p_0)
fitted_0_rho2 = vcat(fitted_0ii_rho2[:,1:7], fitted_0ij_rho2)

### rho3
fitted_ii_0_rho3 = CSV.read("../results/αii_fitted_rho3.csv", DataFrame, header=false)
fitted_ij_0_rho3 = CSV.read("../results/αij_fitted_rho3.csv", DataFrame, header=false)
# fitted = vcat(fitted_ii, fitted_ij)
df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_ii_0_rho3 = DataFrame(fitted_ii_0_rho3, df_names); 
fitted_ij_0_rho3 = DataFrame(fitted_ij_0_rho3, df_names); 

@load "../data/1com0_rho3.jld2" all_ℵii all_ℵij # all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
all_ℵii_0_rho3 = all_ℵii; all_ℵij_0_rho3 = all_ℵij
path_ii_rho3 = ["../results/αii_fitted_rho3.csv"]
path_ij_rho3 = ["../results/αij_fitted_rho3.csv"]
fitted_0ii_rho3 = filter_fitted(path_ii_rho3[1], p_0)
fitted_0ij_rho3 = filter_fitted(path_ij_rho3[1], p_0)
fitted_0_rho3 = vcat(fitted_0ii_rho3[:,1:7], fitted_0ij_rho3)

### rho4 
fitted_ii_0_rho4 = CSV.read("../results/αii_fitted_rho4.csv", DataFrame, header=false)
fitted_ij_0_rho4 = CSV.read("../results/αij_fitted_rho4.csv", DataFrame, header=false)
# fitted = vcat(fitted_ii, fitted_ij)
df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_ii_0_rho4 = DataFrame(fitted_ii_0_rho4, df_names); 
fitted_ij_0_rho4 = DataFrame(fitted_ij_0_rho4, df_names); 

@load "../data/1com0_rho4.jld2" all_ℵii all_ℵij # all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
all_ℵii_0_rho4 = all_ℵii; all_ℵij_0_rho4 = all_ℵij
path_ii_rho4 = ["../results/αii_fitted_rho4.csv"]
path_ij_rho4 = ["../results/αij_fitted_rho4.csv"]
fitted_0ii_rho4 = filter_fitted(path_ii_rho4[1], p_0)
fitted_0ij_rho4 = filter_fitted(path_ij_rho4[1], p_0)
fitted_0_rho4 = vcat(fitted_0ii_rho4[:,1:7], fitted_0ij_rho4)

#################### PLOTTING ##############################
# num_temps = 38; Temp_rich = range(0, num_temps-1, length = num_temps); temp = collect(Temp_rich .+273.15)

f = Figure(fontsize = 30, size = (1800, 1900));
p1 = [PolyElement(color = ("#FF9776", 0.6), strokecolor = "#FF9776", strokewidth = 3)]
p2 = [PolyElement(color = ("#E8C99E", 0.5), strokecolor = "#E8C99E", strokewidth = 3)]
p3 = [PolyElement(color = ("#C25E8B", 0.8), strokecolor = "#C25E8B", strokewidth = 3)]
p4 = [PolyElement(color = ("#82AC6D", 0.5), strokecolor = "#82AC6D", strokewidth = 3)]
p5 = [PolyElement(color = ("#C1C6E8", 0.6), strokecolor = "#C1C6E8", strokewidth = 3)]
p6 = [PolyElement(color = ("#5676A5", 0.8), strokecolor = "#5676A5", strokewidth = 3)]

##### rho2: leaching resource
Label(f[1,0], "Leaching Resource", fontsize = 50, rotation = pi/2)
Box(f[1,1], linestyle = :solid, color = :white)
for i in 1:25
    # Random.seed!(i*5)
    Random.seed!(i); n = rand(1:5)
    ax1 = Axis(f[1,1][all_ind[i][1], all_ind[i][2]], ygridvisible = false, xgridvisible = false)
    hidedecorations!(ax1)
    if i <= 5
        nii_0_rho2 = Int.(range(5*(n-1)+1,5*(n-1)+5,5)[i])
        params_ii_0_rho2 = fitted_ii_0_rho2[Int(nii_0_rho2),1:4]
        αii_0_rho2 = [all_ℵii_0_rho2[t][nii_0_rho2] for t in 1:num_temps]
        pred_ii_0_rho2 = abs.(temp_SS(temp, params_ii_0_rho2))

        # text!(0.0, maximum(abs.(αii)), text = "αᵢᵢ", align = (:left, :top),fontsize = 15)
        scatter!(ax1, Temp_rich, abs.(αii_0_rho2), color = "#5676A5", alpha = 0.9)
        lines!(ax1, Temp_rich, pred_ii_0_rho2, color = ("#FA8328", 0.7), linewidth = 5)
    else 
        nij_0_rho2 = Int.(range(20*(n-1)+1,20*(n-1)+20,20)[i-5])
        params_ij_0_rho2 = fitted_ij_0_rho2[Int(nij_0_rho2),1:4]
        αij_0_rho2 = [all_ℵij_0_rho2[t][nij_0_rho2] for t in 1:num_temps]
        pred_ij_0_rho2 = abs.(temp_SS(temp, params_ij_0_rho2))
        # text!(0.0, maximum(abs.(αij)), text = "αᵢⱼ", align = (:left, :top) ,fontsize = 15)
        scatter!(ax1, Temp_rich, abs.(αij_0_rho2), color = "#5676A5", alpha = 0.9)
        lines!(ax1, Temp_rich, pred_ij_0_rho2, color = ("#015845", 0.7), linewidth = 5)
    end 
end 
ax_diag_0 = Axis(f[1,1], xlabel = "Temperature", ylabel = "|α|", xlabelsize = 35, ylabelsize = 35)
hidedecorations!(ax_diag_0, label = false); hidespines!(ax_diag_0)
lines!(ax_diag_0,[0,1],[1,0], color = (:black, 1.0), linewidth = 3, linestyle = :dash)
Label(f[1,1][1,1, TopLeft()], "(a)")

ax2 = Axis(f[1,2], xlabel = "E", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
xlims!(nothing, 7)
density!(ax2, p_0.E[:,1], label = "Eu", color = ("#FF9776", 0.6), strokewidth = 3, strokecolor = "#FF9776")
lines!(ax2, [median(p_0.E[:,1]), median(p_0.E[:,1])], [0, 1.6], linestyle = :dash, color = ("#AD5525", 0.9), linewidth = 5)
lines!(ax2, [median(p_0.E[:,1]), 3.5],[0.9, 0.9], linestyle = :dot, color = ("#AD5525", 0.9), linewidth = 3)
text!(ax2, 3.5, 0.9, text = "$(round(median(p_0.E[:,1]),digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#AD5525")

density!(ax2, p_0.E[:,2], label = "Em", color = ("#E8C99E", 0.5), strokewidth = 3, strokecolor = "#E8C99E")
lines!(ax2, [median(p_0.E[:,2]), median(p_0.E[:,2])], [0, 1.6], linestyle = :dash, color = ("#F8BA17", 0.9), linewidth = 5)
lines!(ax2, [median(p_0.E[:,2]), 3.5],[1.4, 1.4], linestyle = :dot, color = ("#F8BA17", 0.9), linewidth = 3)
text!(ax2, 3.5, 1.4, text = "$(round(median(p_0.E[:,2]),digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#F8BA17")

density!(ax2, fitted_0_rho2.E, label = "Eα", color = ("#C25E8B", 0.8), strokewidth = 3, strokecolor = "#C25E8B")
lines!(ax2, [median(fitted_0_rho2.E), median(fitted_0_rho2.E)], [0, 1.6], linestyle = :dash, color = ("#601210", 1.0), linewidth = 5)
lines!(ax2, [median(fitted_0_rho2.E), 3.5],[0.6, 0.6], linestyle = :dot, color = ("#601210", 0.9), linewidth = 3)
text!(ax2, 3.5, 0.6, text = "$(round(median(fitted_0_rho2.E) ,digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#601210")

Legend(f[1,2], [p1, p2, p3], tellheight = false, tellwidth = false, [ "Eu", "Em", "Eα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"

Label(f[1,2, TopLeft()], "(b)")

ax3 = Axis(f[1,3], xlabel = "Tₚₖ", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax3, p_0.Tp[:,1] .- 273.15, label = "Tpu", color = ("#82AC6D", 0.5), strokewidth = 3, strokecolor = "#82AC6D")
lines!(ax3, [median(p_0.Tp[:,1])- 273.15, median(p_0.Tp[:,1])- 273.15], [0, 0.13], linestyle = :dash, color = ("#12473D", 0.9), linewidth = 5)
lines!(ax3, [20, median(p_0.Tp[:,1])- 273.15],[0.085, 0.085], linestyle = :dot, color = ("#12473D", 0.9), linewidth = 3)
text!(ax3, 20, 0.085, text = "$(round(median(p_0.Tp[:,1])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#12473D")

density!(ax3, p_0.Tp[:,2] .- 273.15, label = "Tpm", color = ("#C1C6E8", 0.6), strokewidth = 3, strokecolor = "#C1C6E8")
lines!(ax3, [median(p_0.Tp[:,2])- 273.15, median(p_0.Tp[:,2])- 273.15], [0, 0.13], linestyle = :dash, color = ("#9585B4", 0.9), linewidth = 5)
lines!(ax3, [20, median(p_0.Tp[:,2])- 273.15],[0.05, 0.05], linestyle = :dot, color = ("#9585B4", 0.9), linewidth = 3)
text!(ax3, 20, 0.05, text = "$(round(median(p_0.Tp[:,2])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#9585B4")

density!(ax3, fitted_0_rho2.Th .- 273.15, label = "Tpα", color = ("#5676A5", 0.8), strokewidth = 3, strokecolor = "#5676A5")
lines!(ax3, [median(fitted_0_rho2.Th)- 273.15, median(fitted_0_rho2.Th)- 273.15], [0, 0.13], linestyle = :dash, color = ("#0758AE", 0.9), linewidth = 5)
lines!(ax3, [20, median(fitted_0_rho2.Th)- 273.15],[0.095, 0.095], linestyle = :dot, color = ("#0758AE", 0.9), linewidth = 3)
text!(ax3, 20, 0.095, text = "$(round(median(fitted_0_rho2.Th)- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#0758AE")

Legend(f[1,3], [p4, p5, p6], tellheight = false, tellwidth = false, [ "Tpu", "Tpm", "Tpα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,3, TopLeft()], "(c)")

##### rho3: chemostat resource
# f = Figure(fontsize = 30, size = (1800, 2100));

Label(f[2,0], "Chemostat Resource", fontsize = 50, rotation = pi/2)
Box(f[2,1], linestyle = :solid, color = :white)
for i in 1:25
    # Random.seed!(i*5)
    Random.seed!(i); n = rand(1:5)
    ax4 = Axis(f[2,1][all_ind[i][1], all_ind[i][2]], ygridvisible = false, xgridvisible = false)
    hidedecorations!(ax4)
    if i <= 5
        nii_0_rho3 = Int.(range(5*(n-1)+1,5*(n-1)+5,5)[i])
        params_ii_0_rho3 = fitted_ii_0_rho3[Int(nii_0_rho3),1:4]
        αii_0_rho3 = [all_ℵii_0_rho3[t][nii_0_rho3] for t in 1:num_temps]
        pred_ii_0_rho3 = abs.(temp_SS(temp, params_ii_0_rho3))
        # text!(0.0, maximum(abs.(αii)), text = "αᵢᵢ", align = (:left, :top),fontsize = 15)
        scatter!(ax4, Temp_rich, abs.(αii_0_rho3), color = "#5676A5", alpha = 0.9)
        lines!(ax4, Temp_rich, pred_ii_0_rho3, color = ("#FA8328", 0.7), linewidth = 5)
    else 
        nij_0_rho3 = Int.(range(20*(n-1)+1,20*(n-1)+20,20)[i-5])
        params_ij_0_rho3 = fitted_ij_0_rho3[Int(nij_0_rho3),1:4]
        αij_0_rho3 = [all_ℵij_0_rho3[t][nij_0_rho3] for t in 1:num_temps]
        pred_ij_0_rho3 = abs.(temp_SS(temp, params_ij_0_rho3))
        # text!(0.0, maximum(abs.(αij)), text = "αᵢⱼ", align = (:left, :top) ,fontsize = 15)
        scatter!(ax4, Temp_rich, abs.(αij_0_rho3), color = "#5676A5", alpha = 0.9)
        lines!(ax4, Temp_rich, pred_ij_0_rho3, color = ("#015845", 0.7), linewidth = 5)
    end 
end 
ax_diag_0 = Axis(f[2,1], xlabel = "Temperature", ylabel = "|α|", xlabelsize = 35, ylabelsize = 35)
hidedecorations!(ax_diag_0, label = false); hidespines!(ax_diag_0)
lines!(ax_diag_0,[0,1],[1,0], color = (:black, 1.0), linewidth = 3, linestyle = :dash)
Label(f[2,1][1,1, TopLeft()], "(d)")

ax5 = Axis(f[2,2], xlabel = "E", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
xlims!(nothing, 7)
density!(ax5, p_0.E[:,1], label = "Eu", color = ("#FF9776", 0.6), strokewidth = 3, strokecolor = "#FF9776")
lines!(ax5, [median(p_0.E[:,1]), median(p_0.E[:,1])], [0, 1.6], linestyle = :dash, color = ("#AD5525", 0.9), linewidth = 5)
lines!(ax5, [median(p_0.E[:,1]), 3.5],[0.9, 0.9], linestyle = :dot, color = ("#AD5525", 0.9), linewidth = 3)
text!(ax5, 3.5, 0.9, text = "$(round(median(p_0.E[:,1]),digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#AD5525")

density!(ax5, p_0.E[:,2], label = "Em", color = ("#E8C99E", 0.5), strokewidth = 3, strokecolor = "#E8C99E")
lines!(ax5, [median(p_0.E[:,2]), median(p_0.E[:,2])], [0, 1.6], linestyle = :dash, color = ("#F8BA17", 0.9), linewidth = 5)
lines!(ax5, [median(p_0.E[:,2]), 3.5],[1.4, 1.4], linestyle = :dot, color = ("#F8BA17", 0.9), linewidth = 3)
text!(ax5, 3.5, 1.4, text = "$(round(median(p_0.E[:,2]),digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#F8BA17")

density!(ax5, fitted_0_rho3.E, label = "Eα", color = ("#C25E8B", 0.8), strokewidth = 3, strokecolor = "#C25E8B")
lines!(ax5, [median(fitted_0_rho3.E), median(fitted_0_rho3.E)], [0, 1.6], linestyle = :dash, color = ("#601210", 1.0), linewidth = 5)
lines!(ax5, [median(fitted_0_rho3.E), 3.5],[0.6, 0.6], linestyle = :dot, color = ("#601210", 0.9), linewidth = 3)
text!(ax5, 3.5, 0.6, text = "$(round(median(fitted_0_rho3.E) ,digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#601210")

Legend(f[2,2], [p1, p2, p3], tellheight = false, tellwidth = false, [ "Eu", "Em", "Eα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"

Label(f[2,2, TopLeft()], "(e)")

ax6 = Axis(f[2,3], xlabel = "Tₚₖ", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax6, p_0.Tp[:,1] .- 273.15, label = "Tpu", color = ("#82AC6D", 0.5), strokewidth = 3, strokecolor = "#82AC6D")
lines!(ax6, [median(p_0.Tp[:,1])- 273.15, median(p_0.Tp[:,1])- 273.15], [0, 0.14], linestyle = :dash, color = ("#12473D", 0.9), linewidth = 5)
lines!(ax6, [20, median(p_0.Tp[:,1])- 273.15],[0.085, 0.085], linestyle = :dot, color = ("#12473D", 0.9), linewidth = 3)
text!(ax6, 20, 0.085, text = "$(round(median(p_0.Tp[:,1])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#12473D")

density!(ax6, p_0.Tp[:,2] .- 273.15, label = "Tpm", color = ("#C1C6E8", 0.6), strokewidth = 3, strokecolor = "#C1C6E8")
lines!(ax6, [median(p_0.Tp[:,2])- 273.15, median(p_0.Tp[:,2])- 273.15], [0, 0.14], linestyle = :dash, color = ("#9585B4", 0.9), linewidth = 5)
lines!(ax6, [20, median(p_0.Tp[:,2])- 273.15],[0.05, 0.05], linestyle = :dot, color = ("#9585B4", 0.9), linewidth = 3)
text!(ax6, 20, 0.05, text = "$(round(median(p_0.Tp[:,2])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#9585B4")

density!(ax6, fitted_0_rho3.Th .- 273.15, label = "Tpα", color = ("#5676A5", 0.8), strokewidth = 3, strokecolor = "#5676A5")
lines!(ax6, [median(fitted_0_rho3.Th)- 273.15, median(fitted_0_rho3.Th)- 273.15], [0, 0.14], linestyle = :dash, color = ("#0758AE", 0.9), linewidth = 5)
lines!(ax6, [20, median(fitted_0_rho3.Th)- 273.15],[0.095, 0.095], linestyle = :dot, color = ("#0758AE", 0.9), linewidth = 3)
text!(ax6, 20, 0.095, text = "$(round(median(fitted_0_rho3.Th)- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#0758AE")

Legend(f[2,3], [p4, p5, p6], tellheight = false, tellwidth = false, [ "Tpu", "Tpm", "Tpα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[2,3, TopLeft()], "(f)")

################  rho4: self-renewing resource ##### 
# num_temps = 37; Temp_rich = range(0, num_temps-1, length = num_temps); temp = collect(Temp_rich .+273.15)

# f = Figure(fontsize = 30, size = (1800, 1900));

Label(f[3,0], "Self-renewing Resource", fontsize = 50, rotation = pi/2)
Box(f[3,1], linestyle = :solid, color = :white)
for i in 1:25
    # Random.seed!(i*5)
    Random.seed!(i); n = rand(1:5)
    ax7 = Axis(f[3,1][all_ind[i][1], all_ind[i][2]], ygridvisible = false, xgridvisible = false)
    hidedecorations!(ax7)
    if i <= 5
        nii_0_rho4 = Int.(range(5*(n-1)+1,5*(n-1)+5,5)[i])
        params_ii_0_rho4 = fitted_ii_0_rho4[Int(nii_0_rho4),1:4]
        αii_0_rho4 = [all_ℵii_0_rho4[t][nii_0_rho4] for t in 1:num_temps]
        pred_ii_0_rho4 = abs.(temp_SS(temp, params_ii_0_rho4))
        # text!(0.0, maximum(abs.(αii)), text = "αᵢᵢ", align = (:left, :top),fontsize = 15)
        scatter!(ax7, Temp_rich, abs.(αii_0_rho4), color = "#5676A5", alpha = 0.9)
        lines!(ax7, Temp_rich, pred_ii_0_rho4, color = ("#FA8328", 0.7), linewidth = 5)
    else 
        nij_0_rho4 = Int.(range(20*(n-1)+1,20*(n-1)+20,20)[i-5])
        params_ij_0_rho4 = fitted_ij_0_rho4[Int(nij_0_rho4),1:4]
        αij_0_rho4 = [all_ℵij_0_rho4[t][nij_0_rho4] for t in 1:num_temps]
        pred_ij_0_rho4 = abs.(temp_SS(temp, params_ij_0_rho4))
        # text!(0.0, maximum(abs.(αij)), text = "αᵢⱼ", align = (:left, :top) ,fontsize = 15)
        scatter!(ax7, Temp_rich, abs.(αij_0_rho4), color = "#5676A5", alpha = 0.9)
        lines!(ax7, Temp_rich, pred_ij_0_rho4, color = ("#015845", 0.7), linewidth = 5)
    end 
end 
ax_diag_0 = Axis(f[3,1], xlabel = "Temperature", ylabel = "|α|", xlabelsize = 35, ylabelsize = 35)
hidedecorations!(ax_diag_0, label = false); hidespines!(ax_diag_0)
lines!(ax_diag_0,[0,1],[1,0], color = (:black, 1.0), linewidth = 3, linestyle = :dash)
Label(f[3,1][1,1, TopLeft()], "(g)")

ax8 = Axis(f[3,2], xlabel = "E", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
xlims!(nothing, 7)
density!(ax8, p_0.E[:,1], label = "Eu", color = ("#FF9776", 0.6), strokewidth = 3, strokecolor = "#FF9776")
lines!(ax8, [median(p_0.E[:,1]), median(p_0.E[:,1])], [0, 1.6], linestyle = :dash, color = ("#AD5525", 0.9), linewidth = 5)
lines!(ax8, [median(p_0.E[:,1]), 3.5],[0.9, 0.9], linestyle = :dot, color = ("#AD5525", 0.9), linewidth = 3)
text!(ax8, 3.5, 0.9, text = "$(round(median(p_0.E[:,1]),digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#AD5525")

density!(ax8, p_0.E[:,2], label = "Em", color = ("#E8C99E", 0.5), strokewidth = 3, strokecolor = "#E8C99E")
lines!(ax8, [median(p_0.E[:,2]), median(p_0.E[:,2])], [0, 1.6], linestyle = :dash, color = ("#F8BA17", 0.9), linewidth = 5)
lines!(ax8, [median(p_0.E[:,2]), 3.5],[1.4, 1.4], linestyle = :dot, color = ("#F8BA17", 0.9), linewidth = 3)
text!(ax8, 3.5, 1.4, text = "$(round(median(p_0.E[:,2]),digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#F8BA17")

density!(ax8, fitted_0_rho4.E, label = "Eα", color = ("#C25E8B", 0.8), strokewidth = 3, strokecolor = "#C25E8B")
lines!(ax8, [median(fitted_0_rho4.E), median(fitted_0_rho4.E)], [0, 1.6], linestyle = :dash, color = ("#601210", 1.0), linewidth = 5)
lines!(ax8, [median(fitted_0_rho4.E), 3.5],[0.6, 0.6], linestyle = :dot, color = ("#601210", 0.9), linewidth = 3)
text!(ax8, 3.5, 0.6, text = "$(round(median(fitted_0_rho4.E) ,digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#601210")

Legend(f[3,2], [p1, p2, p3], tellheight = false, tellwidth = false, [ "Eu", "Em", "Eα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"

# axislegend(position = :rt)
Label(f[3,2, TopLeft()], "(h)")

ax9 = Axis(f[3,3], xlabel = "Tₚₖ", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax9, p_0.Tp[:,1] .- 273.15, label = "Tpu", color = ("#82AC6D", 0.5), strokewidth = 3, strokecolor = "#82AC6D")
lines!(ax9, [median(p_0.Tp[:,1])- 273.15, median(p_0.Tp[:,1])- 273.15], [0, 1.05], linestyle = :dash, color = ("#12473D", 0.9), linewidth = 5)
lines!(ax9, [20, median(p_0.Tp[:,1])- 273.15],[0.5, 0.5], linestyle = :dot, color = ("#12473D", 0.9), linewidth = 3)
text!(ax9, 20, 0.5, text = "$(round(median(p_0.Tp[:,1])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#12473D")

density!(ax9, p_0.Tp[:,2] .- 273.15, label = "Tpm", color = ("#C1C6E8", 0.6), strokewidth = 3, strokecolor = "#C1C6E8")
lines!(ax9, [median(p_0.Tp[:,2])- 273.15, median(p_0.Tp[:,2])- 273.15], [0, 1.05], linestyle = :dash, color = ("#9585B4", 0.9), linewidth = 5)
lines!(ax9, [20, median(p_0.Tp[:,2])- 273.15],[0.3, 0.3], linestyle = :dot, color = ("#9585B4", 0.9), linewidth = 3)
text!(ax9, 20, 0.3, text = "$(round(median(p_0.Tp[:,2])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#9585B4")

density!(ax9, fitted_0_rho4.Th .- 273.15, label = "Tpα", color = ("#5676A5", 0.8), strokewidth = 3, strokecolor = "#5676A5")
lines!(ax9, [median(fitted_0_rho4.Th)- 273.15, median(fitted_0_rho4.Th)- 273.15], [0, 1.05], linestyle = :dash, color = ("#0758AE", 0.9), linewidth = 5)
lines!(ax9, [20, median(fitted_0_rho4.Th)- 273.15],[0.7, 0.7], linestyle = :dot, color = ("#0758AE", 0.9), linewidth = 3)
text!(ax9, 20, 0.7, text = "$(round(median(fitted_0_rho4.Th)- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#0758AE")

Legend(f[3,3], [p4, p5, p6], tellheight = false, tellwidth = false, [ "Tpu", "Tpm", "Tpα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[3,3, TopLeft()], "(i)")

f

save("../results/TPCα_rho.pdf", f) 
