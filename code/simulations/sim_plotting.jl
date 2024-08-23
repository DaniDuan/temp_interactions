include("./sim_frame.jl");

Eff_results = CSV.read("../results/Eff_results_p0.csv", DataFrame, header=false)
col_names_EF = ["αii", "αii_err", "αij", "αij_err", "αij_d", "αij_d_err", "αij_upper", "αij_upper_err","αij_lower", "αij_lower_err",
                "αii_sur", "αii_sur_err", "αij_sur", "αij_sur_err", "αij_d_sur", "αij_d_sur_err", "αij_upper_sur", "αij_upper_sur_err","αij_lower_sur", "αij_lower_sur_err",
                "r", "r_err", "u", "u_err","m", "m_err", 
                "RO", "RO_err", "ulO", "ulO_err", "estα", "estα_err",
                "RO_sur", "RO_sur_err", "ulO_sur", "ulO_sur_err", "estα_sur", "estα_sur_err",
                "Eu", "Eu_err", "Em", "Em_err", "Eu_sur", "Eu_sur_err", "Em_sur", "Em_sur_err",
                "Tpu", "Tpu_err", "Tpm", "Tpm_err", "Tpu_sur", "Tpu_sur_err", "Tpm_sur", "Tpm_sur_err",
                "eigen", "eigen_err", "stability",
                "diag_dom", "diag_dom_err",
                "Jac_diag", "Jac_diag_err", "radius", "radius_err", "sum_αij", "sum_αij_err"];
rename!(Eff_results, col_names_EF)

N=100
M=50
### Temp params 
# ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
num_temps = 31

### Plots setting ###
Temp_rich = range(0, num_temps-1, length = num_temps)
CairoMakie.activate!(type = "png")
# k = 0.0000862 # Boltzman constant
# x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
# plot(x, log.(abs.(Eff_results.αii)))

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
# ax2 = Axis(f[1,1], ylabel = "log(|αij/αii|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αii, color = ("#FA8328",0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, Eff_results.αii .- Eff_results.αii_err, Eff_results.αii .+ Eff_results.αii_err, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, Eff_results.αij, color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich,  Eff_results.αij .- Eff_results.αij_err, Eff_results.αij .+ Eff_results.αij_err, color = ("#015845", 0.2))
axislegend(position = :rt)
f
# save("../results/a0_T.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|α|)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
# ax2 = Axis(f[1,1], ylabel = "log(|αij/αii|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, log.(abs.(Eff_results.αii)), color = ("#FA8328",0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, log.(abs.(Eff_results.αii .- Eff_results.αii_err)), log.(abs.(Eff_results.αii .+ Eff_results.αii_err)), color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, log.(abs.(Eff_results.αij)), color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich,  log.(abs.(Eff_results.αij .- Eff_results.αij_err)), log.(abs.(Eff_results.αij .+ Eff_results.αij_err)), color = ("#015845", 0.2))
axislegend(position = :lt)
f
# save("../results/a_-1.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50)
lines!(ax, Temp_rich, log.(abs.(Eff_results.αij_d)), color = ("#285C93",1), linewidth = 5, label = "")
band!(ax, Temp_rich, log.(abs.(Eff_results.αij_d .- Eff_results.αij_d_err)) , log.(abs.(Eff_results.αij_d .+ Eff_results.αij_d_err)) , color = ("#285C93", 0.2))
# axislegend(position = :rb)
f
# save("../results/aiiaij_-1.png", f) 

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
# save("../results/CR_RO_-1-1.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Pairwise Interaction", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αij_d, color = ("#285C93",1), linewidth = 5, label = "")
band!(ax1, Temp_rich, Eff_results.αij_d .- Eff_results.αij_d_err , Eff_results.αij_d .+ Eff_results.αij_d_err, color = ("#285C93", 0.2))
lines!(ax2, Temp_rich, Eff_results.estα_sur, color = ("#E17542",1), linewidth = 5, label = "")
band!(ax2, Temp_rich, Eff_results.estα_sur .- Eff_results.estα_sur_err , Eff_results.estα_sur .+ Eff_results.estα_sur_err , color = ("#E17542", 0.2))
l1 = [LineElement(color = ("#285C93",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#E17542", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αij/αii", "ƒc-ƒo"], halign = :center, valign = :top)
f
save("../results/CR_α0.png", f) 

# f = Figure(fontsize = 35, resolution = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Interaction Strength", xlabelsize = 50, ylabelsize = 50)
# lines!(ax, Temp_rich, Eff_results.estα, color = ("#EF8F8C",1), linewidth = 5, label = "")
# band!(ax, Temp_rich, Eff_results.estα .- Eff_results.estα_err , Eff_results.estα .+ Eff_results.estα_err , color = ("#EF8F8C", 0.2))
# # axislegend(position = :rb)
# f
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
lines!(ax2, Temp_rich, abs.(Eff_results.αii), color = ("#EF8F8C", 0.8), linewidth = 5, label = " ")
band!(ax2, Temp_rich,  abs.(Eff_results.αii) .- Eff_results.αii_err, abs.(Eff_results.αii) .+ Eff_results.αii_err, color = ("#EF8F8C", 0.2))
linkxaxes!(ax1,ax2)
# l1 = [LineElement(color = ("#FA8328",0.8), linestyle = nothing, linewidth = 5)]
# l2 = [LineElement(color = ("#015845", 0.8), linestyle = nothing, linewidth = 5)]
# Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
# Label(f[1,1, TopLeft()], "(a)")
f

