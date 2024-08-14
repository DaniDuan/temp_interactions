include("./sim_frame.jl");

# Eff_results = CSV.read("../data/Eff_results0.csv", DataFrame, header=false)
# rename!(Eff_results, col_names_EF)

### Plots setting ###
Temp_rich = range(0, num_temps-1, length = num_temps)
CairoMakie.activate!(type = "png")
# k = 0.0000862 # Boltzman constant
# x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
# plot(x, log.(abs.(Eff_results.αii)))

# data = DataFrame(y = log.(abs.(Eff_results.αii)), x = x);
# Eα = coef(lm(@formula(y ~ x), data))[2]

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|α|)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
# ax2 = Axis(f[1,1], ylabel = "log(|αij/αii|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, abs.(Eff_results.αii), color = ("#FA8328",0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, abs.(Eff_results.αii .- Eff_results.αii_err), abs.(Eff_results.αii .+ Eff_results.αii_err), color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, abs.(Eff_results.αij), color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich,  abs.(Eff_results.αij .- Eff_results.αij_err), abs.(Eff_results.αij .+ Eff_results.αij_err), color = ("#015845", 0.2))
axislegend(position = :lt)
f


f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|α|)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
# ax2 = Axis(f[1,1], ylabel = "log(|αij/αii|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, log.(abs.(Eff_results.αii)), color = ("#FA8328",0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, log.(abs.(Eff_results.αii .- Eff_results.αii_err)), log.(abs.(Eff_results.αii .+ Eff_results.αii_err)), color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, log.(abs.(Eff_results.αij)), color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich,  log.(abs.(Eff_results.αij .- Eff_results.αij_err)), log.(abs.(Eff_results.αij .+ Eff_results.αij_err)), color = ("#015845", 0.2))
axislegend(position = :lt)
f
save("../results/a_-1.png", f) 


f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50)
lines!(ax, Temp_rich, log.(abs.(Eff_results.αij_d)), color = ("#285C93",1), linewidth = 5, label = "")
band!(ax, Temp_rich, log.(abs.(Eff_results.αij_d .- Eff_results.αij_d_err)) , log.(abs.(Eff_results.αij_d .+ Eff_results.αij_d_err)) , color = ("#285C93", 0.2))
# axislegend(position = :rb)
f
save("../results/aiiaij_-1.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Overlap", xlabelsize = 50, ylabelsize = 50)
ax2 = Axis(f[1,1], ylabel = "log(|α|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.ulO, color = ("#FA8328",1), linewidth = 5, label = "Cross-feeding")
band!(ax1, Temp_rich, Eff_results.ulO .- Eff_results.ulO_err , Eff_results.ulO .+ Eff_results.ulO_err , color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, Eff_results.RO, color = ("#015845", 0.6), linewidth = 5, label = "Resource Overlap")
band!(ax1, Temp_rich,  Eff_results.RO .- Eff_results.RO_err, Eff_results.RO .+ Eff_results.RO_err, color = ("#015845", 0.2))
lines!(ax2, Temp_rich, log.(abs.(Eff_results.αii)), color = ("#EF8F8C",0.8), linewidth = 5, label = "αii")
band!(ax2, Temp_rich, log.(abs.(Eff_results.αii .- Eff_results.αii_err)), log.(abs.(Eff_results.αii .+ Eff_results.αii_err)), color = ("#FA8328", 0.2))
lines!(ax2, Temp_rich, log.(abs.(Eff_results.αij)), color = ("#285C93", 0.8), linewidth = 5, label = "αij")
band!(ax2, Temp_rich,  log.(abs.(Eff_results.αij .- Eff_results.αij_err)), log.(abs.(Eff_results.αij .+ Eff_results.αij_err)), color = ("#015845", 0.2))
l1 = [LineElement(color = ("#FA8328",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#015845", 0.8), linestyle = nothing, linewidth = 5)]
l3 = [LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)]
l4 = [LineElement(color = ("#285C93", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2, l3, l4], tellheight = false, tellwidth = false, ["Cross-feeding", "Resource Overlap", "αii" ,"αij"], halign = :left, valign = :top)
f
save("../results/CR_RO_-1-1.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Pairwise Interaction", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αij_d, color = ("#285C93",1), linewidth = 5, label = "")
band!(ax1, Temp_rich, Eff_results.αij_d .- Eff_results.αij_d_err , Eff_results.αij_d .+ Eff_results.αij_d_err, color = ("#285C93", 0.2))
lines!(ax2, Temp_rich, Eff_results.estα, color = ("#E17542",1), linewidth = 5, label = "")
band!(ax2, Temp_rich, Eff_results.estα .- Eff_results.estα_err , Eff_results.estα .+ Eff_results.estα_err , color = ("#E17542", 0.2))
l1 = [LineElement(color = ("#285C93",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#E17542", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αij/αii", "C-R"], halign = :left, valign = :top)
f
save("../results/CR_α0.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Interaction Strength", xlabelsize = 50, ylabelsize = 50)
lines!(ax, Temp_rich, Eff_results.estα, color = ("#EF8F8C",1), linewidth = 5, label = "")
band!(ax, Temp_rich, Eff_results.estα .- Eff_results.estα_err , Eff_results.estα .+ Eff_results.estα_err , color = ("#EF8F8C", 0.2))
# axislegend(position = :rb)
f
# save("../results/IStrength_-1-1.png", f) 

# f = Figure(fontsize = 35, resolution = (1200, 900));
# ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Jacobian", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
# ax2 = Axis(f[1,1], ylabel = "Stability", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false)
# lines!(ax1, Temp_rich, abs.(Eff_results.Jac_diag), color = ("#FA8328",1), linewidth = 5, label = "Diagonal")
# band!(ax1, Temp_rich, abs.(Eff_results.Jac_diag .- Eff_results.Jac_diag_err),  abs.(Eff_results.Jac_diag .+ Eff_results.Jac_diag_err), color = ("#FA8328", 0.2))
# lines!(ax1, Temp_rich, Eff_results.radius, color = ("#015845", 0.6), linewidth = 5, label = "Radius")
# band!(ax1, Temp_rich,  Eff_results.radius .- Eff_results.radius_err, Eff_results.radius .+ Eff_results.radius_err, color = ("#015845", 0.2))
# lines!(ax2, Temp_rich, Eff_results.stability, color = ("#285C93",1), linewidth = 5, label = "Stability")
# # linkxaxes!(ax1,ax2)
# l1 = [LineElement(color = ("#FA8328",0.8), linestyle = nothing, linewidth = 5)]
# l2 = [LineElement(color = ("#015845", 0.8), linestyle = nothing, linewidth = 5)]
# l3 = [LineElement(color = ("#285C93", 0.8), linestyle = nothing, linewidth = 5)]
# Legend(f[1,1], [l1, l2, l3], tellheight = false, tellwidth = false, ["Center", "Radius", "Stability"], halign = :left, valign = :top)
# f
# # save("../results/Jacobian-1.png", f) 

# f = Figure(fontsize = 35, resolution = (1200, 900));
# ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Ratio", xlabelsize = 50, ylabelsize = 50)
# lines!(ax1, Temp_rich, Eff_results.diag_dom, color = ("#FA8328",1), linewidth = 5, label = "Diagonal Dominance")
# band!(ax1, Temp_rich, Eff_results.diag_dom .- Eff_results.diag_dom_err , Eff_results.diag_dom .+ Eff_results.diag_dom_err , color = ("#FA8328", 0.2))
# lines!(ax1, Temp_rich, Eff_results.stability, color = ("#015845",1), linewidth = 5, label = "Stability")
# axislegend(position = :rb)
# f
# # save("../results/Jac_DD_-1.png", f) 

Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = " ", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,1], ylabel = " ", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.u, color = ("#FA8328",0.8), linewidth = 5, label = " ")
band!(ax1, Temp_rich, Eff_results.u .- Eff_results.u_err, Eff_results.u .+ Eff_results.u_err, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, Eff_results.m, color = ("#015845",0.8), linewidth = 5, label = " ")
band!(ax1, Temp_rich, Eff_results.m .- Eff_results.m_err, Eff_results.m .+ Eff_results.m_err, color = ("#015845", 0.2))
lines!(ax2, Temp_rich, Eff_results.α, color = ("#EF8F8C", 0.8), linewidth = 5, label = " ")
band!(ax2, Temp_rich,  Eff_results.α .- Eff_results.α_err, Eff_results.α .+ Eff_results.α_err, color = ("#EF8F8C", 0.2))
linkxaxes!(ax1,ax2)
# l1 = [LineElement(color = ("#FA8328",0.8), linestyle = nothing, linewidth = 5)]
# l2 = [LineElement(color = ("#015845", 0.8), linestyle = nothing, linewidth = 5)]
# Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
# Label(f[1,1, TopLeft()], "(a)")
f
