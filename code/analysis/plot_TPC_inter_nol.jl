include("./fitting.jl");
using ProgressMeter, RCall

N=100
M=50
L = fill(0.0, N)
### Temp params 
num_temps = 38
Tr=273.15+10; Ed=3.5 
x0 = vcat(fill(0.1, N), fill(1, M)) 

condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)

progress = Progress(990; desc="Progress running:")
fitted_all = zeros(Float64, 9900, 6)
for index in 1: 990
    next!(progress)
    @load "../results/inter_ij_0/fit_allij_0_$(index).jld2" fitted
    i = 10*(index-1) + 10
    fitted_all[10*(index-1)+1:i,:] = fitted
end 

df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_all = DataFrame(fitted_all, df_names);
CSV.write("../results/αij_fitted_all0_nol.csv", fitted_all, writeheader=false)

fitted_all_filter = fitted_all[fitted_all.B0 .> 0, :]
fitted_all_filter = fitted_all_filter[fitted_all_filter.E .> 0, :]

##########################################
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
    fitted = rename(fitted, df_names);
    if nrow(fitted) < 101
        fitted.id = tempe_ii
        fitted[!,:"Bu"] = p.B[:,1]
        fitted[!,:"Eu"] = p.E[:,1]
        fitted[!,:"Thu"] = p.Tp[:,1]
    else fitted.id = tempe_ij
    end
    fitted = fitted[fitted.E .> eps(), :]
    fitted = fitted[fitted.B0 .> eps(), :]
    fitted = fitted[fitted.Th .< 323.15, :]
    fitted = fitted[fitted.Th .> 273.15, :]
    fitted = fitted[fitted.r2 .> 0.90, :]
    return fitted
end

idx = collect(CartesianIndices(zeros(Float64, 5, 5)))
ind_diag = diag(idx)
ind_off = [idx[i,j] for i in 1:5 for j in 1:5 if i != j]
all_ind = vcat(ind_diag, ind_off)

########################################
p_0 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=273.15 + 25, ρ_t=[0.0000 0.0000], Tr=Tr, Ed=Ed)

@load "../data/1com_0_nol.jld2" all_ℵii all_ℵij # all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
all_ℵii_0 = all_ℵii; all_ℵij_0 = all_ℵij

fitted_ii_0 = CSV.read("../results/αii_fitted_no_l.csv", DataFrame, header=false)
fitted_ij_0 = CSV.read("../results/αij_fitted_all0_nol.csv", DataFrame, header=false)
# fitted = vcat(fitted_ii, fitted_ij)
df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_ii_0 = DataFrame(fitted_ii_0, df_names)
fitted_ij_0 = DataFrame(fitted_ij_0, df_names)

path_ii = ["../results/αii_fitted_no_l.csv"]
path_ij = ["../results/αij_fitted_all0_nol.csv"]
fitted_0ii = filter_fitted(path_ii[1], p_0)
fitted_0ij = filter_fitted(path_ij[1], p_0)
fitted_0 = vcat(fitted_0ii[:,1:7], fitted_0ij)
# fitted_0r = filter_fitted("../results/r_fitted0.csv", p_0)

f = Figure(fontsize = 30, size = (1800, 700));
#####  ρ = 0 ##### 
# Label(f[1,0], "Minimal Trade-off", fontsize = 50, rotation = pi/2)
Box(f[1,1], linestyle = :solid, color = :white)
for i in 1:25
    # Random.seed!(i*5)
    Random.seed!(i-1); n = rand(1:5)
    ax1 = Axis(f[1,1][all_ind[i][1], all_ind[i][2]], ygridvisible = false, xgridvisible = false)
    hidedecorations!(ax1)
    if i <= 5
        nii_0 = Int.(range(5*(n-1)+1,5*(n-1)+5,5)[i])
        params_ii_0 = fitted_ii_0[Int(nii_0),1:4]
        αii_0 = [all_ℵii_0[t][nii_0] for t in 1:num_temps]
        pred_ii_0 = abs.(temp_SS(temp, params_ii_0))
        # text!(0.0, maximum(abs.(αii)), text = "αᵢᵢ", align = (:left, :top),fontsize = 15)
        scatter!(ax1, Temp_rich, abs.(αii_0), color = "#5676A5", alpha = 0.9)
        lines!(ax1, Temp_rich, pred_ii_0, color = ("#FA8328", 0.7), linewidth = 5)
    else 
        nij_0 = Int.(range(20*(n-1)+1,20*(n-1)+20,20)[i-5])
        params_ij_0 = fitted_ij_0[Int(nij_0),1:4]
        αij_0 = [all_ℵij_0[t][nij_0] for t in 1:num_temps]
        pred_ij_0 = abs.(temp_SS(temp, params_ij_0))
        # text!(0.0, maximum(abs.(αij)), text = "αᵢⱼ", align = (:left, :top) ,fontsize = 15)
        scatter!(ax1, Temp_rich, abs.(αij_0), color = "#5676A5", alpha = 0.9)
        lines!(ax1, Temp_rich, pred_ij_0, color = ("#015845", 0.7), linewidth = 5)
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

density!(ax2, fitted_0.E, label = "Eα", color = ("#C25E8B", 0.8), strokewidth = 3, strokecolor = "#C25E8B")
lines!(ax2, [median(fitted_0.E), median(fitted_0.E)], [0, 1.6], linestyle = :dash, color = ("#601210", 1.0), linewidth = 5)
lines!(ax2, [median(fitted_0.E), 3.5],[0.6, 0.6], linestyle = :dot, color = ("#601210", 0.9), linewidth = 3)
text!(ax2, 3.5, 0.6, text = "$(round(median(fitted_0.E) ,digits = 2)) ev", align = (:left, :center), fontsize = 20, color = "#601210")

p1 = [PolyElement(color = ("#FF9776", 0.6), strokecolor = "#FF9776", strokewidth = 3)]
p2 = [PolyElement(color = ("#E8C99E", 0.5), strokecolor = "#E8C99E", strokewidth = 3)]
p3 = [PolyElement(color = ("#C25E8B", 0.8), strokecolor = "#C25E8B", strokewidth = 3)]
Legend(f[1,2], [p1, p2, p3], tellheight = false, tellwidth = false, [ "Eu", "Em", "Eα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"

# axislegend(position = :rt)
Label(f[1,2, TopLeft()], "(b)")

ax3 = Axis(f[1,3], xlabel = "Tₚₖ", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax3, p_0.Tp[:,1] .- 273.15, label = "Tpu", color = ("#82AC6D", 0.5), strokewidth = 3, strokecolor = "#82AC6D")
lines!(ax3, [median(p_0.Tp[:,1])- 273.15, median(p_0.Tp[:,1])- 273.15], [0, 0.1], linestyle = :dash, color = ("#12473D", 0.9), linewidth = 5)
lines!(ax3, [20, median(p_0.Tp[:,1])- 273.15],[0.095, 0.095], linestyle = :dot, color = ("#12473D", 0.9), linewidth = 3)
text!(ax3, 20, 0.095, text = "$(round(median(p_0.Tp[:,1])- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#12473D")

density!(ax3, p_0.Tp[:,2] .- 273.15, label = "Tpm", color = ("#C1C6E8", 0.6), strokewidth = 3, strokecolor = "#C1C6E8")
lines!(ax3, [median(p_0.Tp[:,2])- 273.15, median(p_0.Tp[:,2])- 273.15], [0, 0.1], linestyle = :dash, color = ("#9585B4", 0.9), linewidth = 5)
lines!(ax3, [median(p_0.Tp[:,2])- 273.15, 48],[0.05, 0.05], linestyle = :dot, color = ("#9585B4", 0.9), linewidth = 3)
text!(ax3, 48, 0.05, text = "$(round(median(p_0.Tp[:,2])- 273.15 ,digits = 2)) °C", align = (:left, :center), fontsize = 20, color = "#9585B4")

density!(ax3, fitted_0.Th .- 273.15, label = "Tpα", color = ("#5676A5", 0.8), strokewidth = 3, strokecolor = "#5676A5")
lines!(ax3, [median(fitted_0.Th)- 273.15, median(fitted_0.Th)- 273.15], [0, 0.1], linestyle = :dash, color = ("#0758AE", 0.9), linewidth = 5)
lines!(ax3, [20, median(fitted_0.Th)- 273.15],[0.085, 0.085], linestyle = :dot, color = ("#0758AE", 0.9), linewidth = 3)
text!(ax3, 20, 0.085, text = "$(round(median(fitted_0.Th)- 273.15 ,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#0758AE")

p4 = [PolyElement(color = ("#82AC6D", 0.5), strokecolor = "#82AC6D", strokewidth = 3)]
p5 = [PolyElement(color = ("#C1C6E8", 0.6), strokecolor = "#C1C6E8", strokewidth = 3)]
p6 = [PolyElement(color = ("#5676A5", 0.8), strokecolor = "#5676A5", strokewidth = 3)]
Legend(f[1,3], [p4, p5, p6], tellheight = false, tellwidth = false, [ "Tpu", "Tpm", "Tpα"], halign = :right, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,3, TopLeft()], "(c)")

f
save("../results/TPCα_nol.pdf", f) 
