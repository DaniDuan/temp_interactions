include("./sim_frame.jl")
using ProgressMeter, RCall, Glob, ColorSchemes
num_temps = 31
N=100; M=50
Temp_rich = range(0, num_temps-1, length = num_temps)

############## collecting results ##############
path_0 = glob("Eff_iters*", "../data/Eff_p0_Umatrix/")
progress = Progress(length(path_0)*num_temps; desc="Progress running:")
all_ii_collect_0 = Vector{Vector{Float64}}(); all_ij_collect_0 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    α_ij =  Float64[]; α_ii = Float64[];
    for i in 1:length(path_0)
        @load path_0[i]  all_ℵii all_ℵij #all_sur rich
        append!(α_ii, all_ℵii[j]); append!(α_ij, all_ℵij[j]); 
        next!(progress)
    end 
    push!(all_ii_collect_0, α_ii); push!(all_ij_collect_0, α_ij);
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

path_1 = glob("Eff_iters*", "../data/Eff_p-1_Umatrix/")
progress = Progress(length(path_1)*num_temps; desc="Progress running:")
all_ii_collect_1 = Vector{Vector{Float64}}(); all_ij_collect_1 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    α_ij =  Float64[]; α_ii = Float64[];
    for i in 1:length(path_1)
        @load path_1[i]  all_ℵii all_ℵij #all_sur rich
        append!(α_ii, all_ℵii[j]); append!(α_ij, all_ℵij[j]); 
        next!(progress)
    end 
    push!(all_ii_collect_1, α_ii); push!(all_ij_collect_1, α_ij);
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

############## analysing ##############
αii_0 = [mean(all_ii_collect_0[t]) for t in 1: num_temps]
αii_err_0 = [std(all_ii_collect_0[t])/sqrt(length(all_ii_collect_0[t])) for t in 1: num_temps]
αij_0 = [mean(all_ij_collect_0[t]) for t in 1: num_temps]
αij_err_0 = [std(all_ij_collect_0[t])/sqrt(length(all_ij_collect_0[t])) for t in 1: num_temps]

αii_1 = [mean(all_ii_collect_1[t]) for t in 1: num_temps]
αii_err_1 = [std(all_ii_collect_1[t])/sqrt(length(all_ii_collect_1[t])) for t in 1: num_temps]
αij_1 = [mean(all_ij_collect_1[t]) for t in 1: num_temps]
αij_err_1 = [std(all_ij_collect_1[t])/sqrt(length(all_ij_collect_1[t])) for t in 1: num_temps]

############## plotting ##############
cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])
f = Figure(fontsize = 30, size = (1800, 1200));
Label(f[1,0], "Minimal Trade-off", fontsize = 50, rotation = pi/2)
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, αii_0, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax1, Temp_rich, αii_0 .- αii_err_0, αii_0 .+ αii_err_0, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, αij_0, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
band!(ax1, Temp_rich,  αij_0 .- αij_err_0, αij_0 .+ αij_err_0, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[1,1, TopLeft()], "(a)")
ax2 = Axis(f[1,2], limits = ((-22.0, 3.0), nothing), xlabel = L"log(|α_{ii}|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
for i in 1: 31
    density!(ax2, log.(abs.(all_ii_collect_0[i])), color = (cs[i], 0.2), strokewidth = 3, strokecolor = (cs[i], 0.7))
end 
Label(f[1,2, TopLeft()], "(b)")
ax3 = Axis(f[1,3], limits = ((-22.0, 3.0), nothing), xlabel = L"log(|α_{i≠j}|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
for i in 1: 31
    density!(ax3, log.(abs.(all_ij_collect_0[i])), color = (cs[i], 0.3), strokewidth = 3, strokecolor = (cs[i], 0.5))
end 
Colorbar(f[1,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[1,3, TopLeft()], "(c)")

Label(f[2,0], "Maximal Trade-off", fontsize = 50, rotation = pi/2)
ax4 = Axis(f[2,1],
    xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax4, Temp_rich, αii_1, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax4, Temp_rich, αii_1 .- αii_err_1, αii_1 .+ αii_err_1, color = ("#FA8328", 0.2))
lines!(ax4, Temp_rich, αij_1, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
band!(ax4, Temp_rich,  αij_1 .- αij_err_1, αij_1 .+ αij_err_1, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[2,1, TopLeft()], "(d)")
ax5 = Axis(f[2,2], limits = ((-22.0, 3.0), nothing), xlabel = L"log(|α_{ii}|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
for i in 1: 31
    density!(ax5, log.(abs.(all_ii_collect_1[i])), color = (cs[i], 0.2), strokewidth = 3, strokecolor = (cs[i], 0.7))
end 
Label(f[2,2, TopLeft()], "(e)")
ax6 = Axis(f[2,3], limits = ((-22.0, 3.0), nothing), xlabel = L"log(|α_{i≠j}|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
for i in 1: 31
    density!(ax6, log.(abs.(all_ij_collect_1[i])), color = (cs[i], 0.3), strokewidth = 3, strokecolor = (cs[i], 0.5))
end 
Colorbar(f[2,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[2,3, TopLeft()], "(f)")

f

save("../results/distα.pdf", f) 

# repeating_ii = [vcat([repeat([all_ii_collect[t][i]], 99) for i in 1:length(all_ii_collect[t])]...) for t in 1:num_temps]
# cor_ij_ii = [cor(repeating_ii[t], all_ij_collect[t]) for t in 1:num_temps]
# plot([mean(all_ij_collect[i]) for i in 1:num_temps])
# plot(cor_ij_ii)
