include("./sim_frame.jl")
using ProgressMeter, RCall, Glob, ColorSchemes, KernelDensity, Colors
num_temps = 31
N=100; M=50
Temp_rich = range(0, num_temps-1, length = num_temps)

@load "../data/summary_temp_inter_rho.jld2" all_rich_collect_1 all_ii_collect_1 all_ij_collect_1 all_rich_collect_2 all_ii_collect_2 all_ij_collect_2 all_rich_collect_3 all_ii_collect_3 all_ij_collect_3 all_rich_collect_4 all_ii_collect_4 all_ij_collect_4
############## analysing ##############
αii_ρ1 = [mean(all_ii_collect_1[t]) for t in 1: num_temps]
αii_err_ρ1 = [std(all_ii_collect_1[t])/sqrt(length(all_ii_collect_1[t])) for t in 1: num_temps]
αij_ρ1 = [mean(all_ij_collect_1[t]) for t in 1: num_temps]
αij_err_ρ1 = [std(all_ij_collect_1[t])/sqrt(length(all_ij_collect_1[t])) for t in 1: num_temps]

αii_ρ2 = [mean(all_ii_collect_2[t]) for t in 1: num_temps]
αii_err_ρ2 = [std(all_ii_collect_2[t])/sqrt(length(all_ii_collect_2[t])) for t in 1: num_temps]
αij_ρ2 = [mean(all_ij_collect_2[t]) for t in 1: num_temps]
αij_err_ρ2 = [std(all_ij_collect_2[t])/sqrt(length(all_ij_collect_2[t])) for t in 1: num_temps]

αii_ρ3 = [mean(all_ii_collect_3[t]) for t in 1: num_temps]
αii_err_ρ3 = [std(all_ii_collect_3[t])/sqrt(length(all_ii_collect_3[t])) for t in 1: num_temps]
αij_ρ3 = [mean(all_ij_collect_3[t]) for t in 1: num_temps]
αij_err_ρ3 = [std(all_ij_collect_3[t])/sqrt(length(all_ij_collect_3[t])) for t in 1: num_temps]

αii_ρ4 = [mean(all_ii_collect_4[t]) for t in 1: num_temps]
αii_err_ρ4 = [std(all_ii_collect_4[t])/sqrt(length(all_ii_collect_4[t])) for t in 1: num_temps]
αij_ρ4 = [mean(all_ij_collect_4[t]) for t in 1: num_temps]
αij_err_ρ4 = [std(all_ij_collect_4[t])/sqrt(length(all_ij_collect_4[t])) for t in 1: num_temps]

############## plotting ##############
cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])

# f = Figure(fontsize = 30, size = (1800, 2600));
f = Figure(fontsize = 30, size = (1800, 2000));

Label(f[1,0], "Constant Resource", fontsize = 50, rotation = pi/2)
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α",
    xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, αii_ρ1, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax1, Temp_rich, αii_ρ1 .- αii_err_ρ1, αii_ρ1 .+ αii_err_ρ1, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, αij_ρ1, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
band!(ax1, Temp_rich,  αij_ρ1 .- αij_err_ρ1, αij_ρ1 .+ αij_err_ρ1, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[1,1, TopLeft()], "(a)")

ax2 = Axis(f[1,2], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{ii}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ii_collect_1[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax2, kd.x, kd.density, color = (cs[i], 0.7), linewidth = 3)
        poly!(ax2, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[1,2, TopLeft()], "(b)")

ax3 = Axis(f[1,3], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{i≠j}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ij_collect_1[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax3, kd.x, kd.density, color = (cs[i], 0.7), linewidth = 3)
        poly!(ax3, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Colorbar(f[1,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[1,3, TopLeft()], "(c)")

########### LEACHING RESOURCE 
Label(f[2,0], "Leaching Resource", fontsize = 50, rotation = pi/2)
ax4 = Axis(f[2,1], xlabel = "Temperature (°C)", ylabel = "α",
    xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax4, Temp_rich, αii_ρ2, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax4, Temp_rich, αii_ρ2 .- αii_err_ρ2, αii_ρ2 .+ αii_err_ρ2, color = ("#FA8328", 0.2))
lines!(ax4, Temp_rich, αij_ρ2, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
band!(ax4, Temp_rich,  αij_ρ2 .- αij_err_ρ2, αij_ρ2 .+ αij_err_ρ2, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[2,1, TopLeft()], "(d)")

ax5 = Axis(f[2,2], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{ii}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ii_collect_2[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax5, kd.x, kd.density, color = (cs[i], 0.7), linewidth = 3)
        poly!(ax5, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[2,2, TopLeft()], "(e)")

ax6 = Axis(f[2,3], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{i≠j}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ij_collect_2[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax6, kd.x, kd.density, color = (cs[i], 0.7), linewidth = 3)
        poly!(ax6, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Colorbar(f[2,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[2,3, TopLeft()], "(f)")

########### CHEMOSTAT RESOURCE 
Label(f[3,0], "Chemostat Resource", fontsize = 50, rotation = pi/2)
ax7 = Axis(f[3,1], xlabel = "Temperature (°C)", ylabel = "α",
    xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax7, Temp_rich, αii_ρ3, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax7, Temp_rich, αii_ρ3 .- αii_err_ρ3, αii_ρ3 .+ αii_err_ρ3, color = ("#FA8328", 0.2))
lines!(ax7, Temp_rich, αij_ρ3, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
band!(ax7, Temp_rich,  αij_ρ3 .- αij_err_ρ3, αij_ρ3 .+ αij_err_ρ3, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[3,1, TopLeft()], "(g)")

ax8 = Axis(f[3,2], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{ii}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ii_collect_3[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax8, kd.x, kd.density, color = (cs[i], 0.7), linewidth = 3)
        poly!(ax8, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[3,2, TopLeft()], "(h)")

ax9 = Axis(f[3,3], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{i≠j}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ij_collect_3[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax9, kd.x, kd.density, color = (cs[i], 0.7), linewidth = 3)
        poly!(ax9, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Colorbar(f[3,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[3,3, TopLeft()], "(i)")

Label(f[0, :], "With Leakage (L = 0.3)", fontsize = 60,
    tellwidth = false, font = :bold, justification = :center)

    # f = Figure(fontsize = 30, size = (1800, 2000));
# Label(f[4,0], "Self-renewing Resource", fontsize = 50, rotation = pi/2)
# ax10 = Axis(f[4,1],
#     xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
# lines!(ax10, Temp_rich[1:31], αii_ρ4[1:31], color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
# band!(ax10, Temp_rich[1:31], αii_ρ4[1:31] .- αii_err_ρ4[1:31], αii_ρ4[1:31] .+ αii_err_ρ4[1:31], color = ("#FA8328", 0.2))
# lines!(ax10, Temp_rich[1:31], αij_ρ4[1:31], color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
# band!(ax10, Temp_rich[1:31],  αij_ρ4[1:31] .- αij_err_ρ4[1:31], αij_ρ4[1:31] .+ αij_err_ρ4[1:31], color = ("#015845", 0.2))
# axislegend(position = :lb)
# Label(f[4,1, TopLeft()], "(j)")
# ax11 = Axis(f[4,2], limits = ((-22.0, 3.0), nothing), xlabel = L"log(|α_{ii}|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
# for i in 1: 31
#     clean_data = log.(abs.(all_ii_collect_4[i]))
#     clean_data = clean_data[isfinite.(clean_data)]
#     density!(ax11, clean_data, color = (cs[i], 0.2), strokewidth = 3, strokecolor = (cs[i], 0.7))
# end 
# Label(f[4,2, TopLeft()], "(k)")
# ax12 = Axis(f[4,3], limits = ((-22.0, 3.0), nothing), xlabel = L"log(|α_{i≠j}|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35)
# for i in 1: 31
#     clean_data = log.(abs.(all_ij_collect_4[i]))
#     clean_data = clean_data[isfinite.(clean_data)]
#     density!(ax12, clean_data, color = (cs[i], 0.3), strokewidth = 3, strokecolor = (cs[i], 0.5))
# end 
# Colorbar(f[4,4], colorrange = [0, 25], colormap = cs[1:31], label = "Temperature")
# Label(f[4,3, TopLeft()], "(l)")

Label(f[0, :], "With Leakage (L = 0.3)", fontsize = 60, tellwidth = false, font = :bold, justification = :center)

f
save("../results/distα_rho_v2.pdf", f) 


#####################################################
cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])

f = Figure(fontsize = 30, size = (1800, 700));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α",
    xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, αii_ρ1, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax1, Temp_rich, αii_ρ1 .- αii_err_ρ1, αii_ρ1 .+ αii_err_ρ1, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, αij_ρ1, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i≠j}")
band!(ax1, Temp_rich, αij_ρ1 .- αij_err_ρ1, αij_ρ1 .+ αij_err_ρ1, color = ("#015845", 0.2))
axislegend(position = :lb)
Label(f[1,1, TopLeft()], "(a)")
ax2 = Axis(f[1,2], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{ii}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ii_collect_1[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax2, kd.x, kd.density, color = (cs[i], 0.6), linewidth = 3)
        poly!(ax2, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[1,2, TopLeft()], "(b)")
ax3 = Axis(f[1,3], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{i≠j}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ij_collect_1[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax3, kd.x, kd.density, color = (cs[i], 0.6), linewidth = 3)
        poly!(ax3, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[1,3, TopLeft()], "(c)")
Colorbar(f[1,4], colorrange = [0, 30], colormap = cs, label = "Temperature")
f
save("../results/distα_rho_v3.pdf", f) 


#####################################################
cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])

f = Figure(fontsize = 30, size = (1800, 700));
ax2 = Axis(f[1,1], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{ii}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ii_collect_1[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax2, kd.x, kd.density, color = (cs[i], 0.6), linewidth = 3)
        poly!(ax2, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[1,1, TopLeft()], "(a)")
ax3 = Axis(f[1,2], limits = ((-22.0, 3.0), nothing),
    xlabel = L"log(|α_{j≠i}|)", ylabel = "Density",
    xlabelsize = 35, ylabelsize = 35)
for i in 1:31
    clean_data = log.(abs.(all_ij_collect_1[i]))
    clean_data = clean_data[isfinite.(clean_data)]
    if !isempty(clean_data)
        kd = kde(clean_data; npoints = 300, boundary = (-22.0, 3.0))
        lines!(ax3, kd.x, kd.density, color = (cs[i], 0.6), linewidth = 3)
        poly!(ax3, kd.x, kd.density, color = (cs[i], 0.2))
    end
end
Label(f[1,2, TopLeft()], "(b)")
Colorbar(f[1,3], colorrange = [0, 30], colormap = cs, label = "Temperature")

ax1 = Axis(f[1,4], xlabel = "Temperature (°C)", ylabel = L"\text{Mean Interaction Strength (} \bar{\alpha} \pm \text{SE} \text{)}",
    xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, αii_ρ1, color = ("#FA8328", 0.8), linewidth = 5, label = L"\bar{\alpha}_{ii}")
band!(ax1, Temp_rich, αii_ρ1 .- αii_err_ρ1, αii_ρ1 .+ αii_err_ρ1, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, αij_ρ1, color = ("#015845", 0.8), linewidth = 5, label = L"\bar{\alpha}_{j≠i}")
band!(ax1, Temp_rich, αij_ρ1 .- αij_err_ρ1, αij_ρ1 .+ αij_err_ρ1, color = ("#015845", 0.2))
axislegend(position = :rt)
Label(f[1,4, TopLeft()], "(c)")

f
save("../results/distα_rho_v4.svg", f) 


##############################################
αij_ρ1_var = [mean([var(log.(abs.(all_ij_collect_1[t][(i-1)*9900+1:i*9900]))) for i in 1:999]) for t in 1: num_temps]
αij_err_ρ1_var = [std([var(log.(abs.(all_ij_collect_1[t][(i-1)*9900+1:i*9900]))) for i in 1:999])/sqrt(999) for t in 1: num_temps]
rich_ρ1 = [mean(all_rich_collect_1[t]) for t in 1: num_temps]
rich_err_ρ1 = [std(all_rich_collect_1[t])/sqrt(length(all_rich_collect_1[t])) for t in 1: num_temps]

f = Figure(fontsize = 35, size = (1000, 800));
# Label(f[1,0], "With Leakage (L = 0.3)", fontsize = 50, rotation = pi/2)
# Label(f[0,1], "Minimal Trade-off", fontsize = 50)
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Variation in α", xlabelsize = 45, ylabelsize = 45, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Richness", yaxisposition = :right, yticklabelalign = (:left, :center), xlabelsize = 45, ylabelsize = 45, ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
# hidedecorations!(ax1, grid = false, ticks = true, ticklabels = true)
lines!(ax1, Temp_rich, αij_ρ1_var, color = ("#FA8328",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, αij_ρ1_var .- αij_err_ρ1_var , αij_ρ1_var.+ αij_err_ρ1_var, color = ("#FA8328", 0.3))
lines!(ax2, Temp_rich, rich_ρ1, color =( "#015845", 0.9), linewidth = 5, label = "")
band!(ax2, Temp_rich, rich_ρ1 .- rich_err_ρ1 , rich_ρ1.+ rich_err_ρ1, color = ("#015845", 0.3))
linkxaxes!(ax1,ax2)
# lines!(ax1, [0, 30], [1, 1], linestyle = :dash, color = ("#4F363E", 0.9), linewidth = 2)
# text!(ax1, 0, 1.05, text = L"↑ α_{j≠i} > α_{ii}", align = (:left, :center),fontsize = 30)
# text!(ax1, 0, 0.95, text = L"↓ α_{j≠i} < α_{ii}", align = (:left, :center),fontsize = 30)
l1 = [LineElement(color = ("#FA8328", 0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#015845", 0.9), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ L"var(log(|α_{j≠i; j = 1,...,N}|))", "Richness"], halign = :left, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/var_rich_rho.pdf", f) 

f = Figure(fontsize = 35, size = (1000, 800));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Richness", xlabelsize = 45, ylabelsize = 45, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Variation in α", yaxisposition = :right, yticklabelalign = (:left, :center), xlabelsize = 45, ylabelsize = 45, ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
lines!(ax1, Temp_rich, rich_ρ1, color =( "#015845", 0.9), linewidth = 5, label = "")
band!(ax1, Temp_rich, rich_ρ1 .- rich_err_ρ1 , rich_ρ1.+ rich_err_ρ1, color = ("#015845", 0.3))
lines!(ax2, Temp_rich, αij_ρ1_var, color = ("#FA8328",0.8), linewidth = 5, label = "")
band!(ax2, Temp_rich, αij_ρ1_var .- αij_err_ρ1_var , αij_ρ1_var.+ αij_err_ρ1_var, color = ("#FA8328", 0.3))
linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#FA8328", 0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#015845", 0.9), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ L"var(log(|α_{j≠i; j = 1,...,N}|))", "Richness"], halign = :left, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/var_rich_rho_v1.pdf", f) 


f = Figure(fontsize = 75, size = (1000, 1000));
ax1 = Axis(f[1,1], xlabel = "T (°C)", ylabel = L"Var(log(α_{j≠i}))", xlabelsize = 100, ylabelsize = 100, ygridvisible = false, xgridvisible = false)
lines!(ax1, Temp_rich, αij_ρ1_var, color = ("#FA8328",0.8), linewidth = 8, label = "")
band!(ax1, Temp_rich, αij_ρ1_var .- αij_err_ρ1_var , αij_ρ1_var.+ αij_err_ρ1_var, color = ("#FA8328", 0.3))
f
save("../results/var_rho.svg", f) 



############################
