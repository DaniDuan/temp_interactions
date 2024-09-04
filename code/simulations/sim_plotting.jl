include("./sim_frame.jl");

Eff_results = CSV.read("../results/Eff_results_p0_new.csv", DataFrame, header=false)[:, 1:65]
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
rename!(Eff_results, col_names_EF)
@load "../results/Feas_CR_dist_p0_new.jld2" all_Rrela_collect all_Crela_collect all_R_collect all_C_collect all_ii_collect all_ij_collect all_ii_sur_collect all_ij_sur_collect all_r_collect all_r_sur_collect

Temp_rich = range(0, num_temps-1, length = num_temps)
temp = collect(Temp_rich .+ 273.15)
temp_R = vcat([repeat([Temp_rich[t]], length(all_R_collect[t])) for t in 1:num_temps]...)
temp_C = vcat([repeat([Temp_rich[t]], length(all_C_collect[t])) for t in 1:num_temps]...)
# allR = vcat(all_R_collect...)
# allC = vcat(all_C_collect...)
meanR = [mean(all_R_collect[t]) for t in 1:num_temps]
R_err = [std(all_R_collect[t])/sqrt(length(all_R_collect[t])) for t in 1:num_temps]
meanC = [mean(all_C_collect[t]) for t in 1:num_temps]
C_err = [std(all_C_collect[t])/sqrt(length(all_C_collect[t])) for t in 1:num_temps]

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Growth Rate", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Biomass", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
lines!(ax1, Temp_rich, Eff_results.r_sur, color = ("#82AC6D", 0.9), linewidth = 5, label = "r_sur")
band!(ax1, Temp_rich, Eff_results.r_sur .- Eff_results.r_sur_err, Eff_results.r_sur .+ Eff_results.r_sur_err, color = ("#82AC6D", 0.5))
lines!(ax2, Temp_rich, meanC, color = ("#0758AE", 0.8), linewidth = 5, label = "biomass")
band!(ax2, Temp_rich, meanC .- C_err, meanC .+ C_err, color = ("#0758AE", 0.5))
l1 = [LineElement(color = ("#82AC6D", 0.9), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#0758AE", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "survivor r", "biomass"], halign = :center, valign = :top)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/biomass_r-1.pdf", f) 

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Resource Abundance", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
# ax3 = Axis(f[1,1], ylabel = "", ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
# hidespines!(ax3)
hidedecorations!(ax3, grid = false, ticks = true, ticklabels = true)
lines!(ax1, Temp_rich, Eff_results.αij_d, color = ("#376298",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, Eff_results.αij_d .- Eff_results.αij_d_err , Eff_results.αij_d.+ Eff_results.αij_d_err, color = ("#376298", 0.3))
lines!(ax2, Temp_rich, meanR, color =( "#F8BA17", 0.9), linewidth = 5, label = "")
band!(ax2, Temp_rich, meanR .- R_err , meanR .+ R_err , color = ("#F8BA17", 0.5))
# lines!(ax3, Temp_rich, Eff_results.estα, color = ("#9A2B1A",0.7), linewidth = 5, label = "")
# band!(ax3, Temp_rich, Eff_results.estα .- Eff_results.estα_err , Eff_results.estα .+ Eff_results.estα_err , color = ("#9A2B1A", 0.5))
linkxaxes!(ax1,ax2, ax3)
l1 = [LineElement(color = ("#376298", 0.8), linestyle = nothing, linewidth = 5)]
# l2 = [LineElement(color = ("#9A2B1A", 0.7), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#F8BA17", 0.9), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αij/αii", "Resource"], halign = :center, valign = :top, framevisible = false) # "ƒc-ƒo"
# Label(f[1,1, TopLeft()], "(a)")
f
save("../results/αR0.pdf", f) 

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Resource Abundance", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Consumer Abundance", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
lines!(ax1, Temp_rich, meanR, color =( "#F8BA17", 0.9), linewidth = 5, label = "")
band!(ax1, Temp_rich, meanR .- R_err , meanR .+ R_err , color = ("#F8BA17", 0.5))
lines!(ax2,Temp_rich, meanC, color = ("#0758AE", 0.7), linewidth = 5, label = "")
band!(ax2, Temp_rich, meanC .- C_err , meanC .+ C_err , color = ("#0758AE", 0.5))
linkxaxes!(ax1,ax2)
p1 = [LineElement(color = ("#F8BA17", 0.9), linestyle = nothing, linewidth = 5)]
p2 = [LineElement(color = ("#0758AE", 0.7), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [p1, p2], tellheight = false, tellwidth = false, ["Resource", "Consumer"], halign = :right, valign = :top)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/CR0.pdf", f) 




################################
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
CairoMakie.activate!(type = "pdf")
# k = 0.0000862 # Boltzman constant
# x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
# plot(x, log.(abs.(Eff_results.αii)))

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
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
f
save("../results/a0_T.pdf", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α_sur", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
# ax2 = Axis(f[1,1], ylabel = "log(|αij/αii|)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αii_sur, color = ("#FA8328",0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, Eff_results.αii_sur .- Eff_results.αii_sur_err, Eff_results.αii_sur .+ Eff_results.αii_sur_err, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich,Eff_results.αij_sur, color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich, Eff_results.αij_sur .- Eff_results.αij_sur_err, Eff_results.αij_sur .+ Eff_results.αij_sur_err, color = ("#015845", 0.2))
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

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Pairwise Interaction", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
lines!(ax1, Temp_rich, Eff_results.αij_d, color = ("#376298",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, Eff_results.αij_d .- Eff_results.αij_d_err , Eff_results.αij_d.+ Eff_results.αij_d_err, color = ("#376298", 0.3))
lines!(ax2, Temp_rich, Eff_results.estα, color = ("#9A2B1A",0.7), linewidth = 5, label = "")
band!(ax2, Temp_rich, Eff_results.estα .- Eff_results.estα_err , Eff_results.estα .+ Eff_results.estα_err , color = ("#9A2B1A", 0.3))
l1 = [LineElement(color = ("#376298", 0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#9A2B1A", 0.7), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αij/αii", "ƒc-ƒo"], halign = :center, valign = :top)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/CR_α0.pdf", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αij/αii", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Pairwise Interaction", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αij_d_sur, color = ("#285C93",1), linewidth = 5, label = "")
band!(ax1, Temp_rich, Eff_results.αij_d_sur .- Eff_results.αij_d_sur_err , Eff_results.αij_d_sur.+ Eff_results.αij_d_sur_err, color = ("#285C93", 0.2))
lines!(ax2, Temp_rich, Eff_results.estα_sur, color = ("#E17542",1), linewidth = 5, label = "")
band!(ax2, Temp_rich, Eff_results.estα_sur .- Eff_results.estα_sur_err , Eff_results.estα_sur .+ Eff_results.estα_sur_err , color = ("#E17542", 0.2))
l1 = [LineElement(color = ("#285C93",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#E17542", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αij/αii", "ƒc-ƒo"], halign = :center, valign = :top)
f

# f = Figure(fontsize = 35, resolution = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Interaction Strength", xlabelsize = 50, ylabelsize = 50)
# lines!(ax, Temp_rich, Eff_results.estα, color = ("#EF8F8C",1), linewidth = 5, label = "")
# band!(ax, Temp_rich, Eff_results.estα .- Eff_results.estα_err , Eff_results.estα .+ Eff_results.estα_err , color = ("#EF8F8C", 0.2))
# # axislegend(position = :rb)
# f
# save("../results/IStrength_-1-1.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Jacobian", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,1], ylabel = "Stability", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, abs.(Eff_results.Jac_diag), color = ("#FA8328",1), linewidth = 5, label = "Diagonal")
band!(ax1, Temp_rich, abs.(Eff_results.Jac_diag .- Eff_results.Jac_diag_err),  abs.(Eff_results.Jac_diag .+ Eff_results.Jac_diag_err), color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, Eff_results.radius, color = ("#015845", 0.6), linewidth = 5, label = "Radius")
band!(ax1, Temp_rich,  Eff_results.radius .- Eff_results.radius_err, Eff_results.radius .+ Eff_results.radius_err, color = ("#015845", 0.2))
lines!(ax2, Temp_rich, Eff_results.stability, color = ("#285C93",1), linewidth = 5, label = "Stability")
# linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#FA8328",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#015845", 0.8), linestyle = nothing, linewidth = 5)]
l3 = [LineElement(color = ("#285C93", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2, l3], tellheight = false, tellwidth = false, ["Center", "Radius", "Stability"], halign = :left, valign = :top)
f
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
f = Figure(fontsize = 35, size = (1200, 900));
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
