include("./fitting.jl");
using ProgressMeter, RCall

N=100
M=50
L = fill(0.3, N)
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
    @load "../results/fit_allij-1_new/fit_allij-1_$(index).jld2" fitted
    i = 10*(index-1) + 10
    fitted_all[10*(index-1)+1:i,:] = fitted
end 

df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_all = DataFrame(fitted_all, df_names);
CSV.write("../results/αij_fitted_all-1_new.csv", fitted_all, writeheader=false)

# fitted_all = CSV.read("../results/αij_fitted_all_0_new.csv", DataFrame, header=false)
# df_names = ["B0","E","Th","Ed","AIC","r2"]
# fitted_all = DataFrame(fitted_all, df_names);

fitted_all_filter = fitted_all[fitted_all.B0 .> 0, :]
fitted_all_filter = fitted_all_filter[fitted_all_filter.E .> 0, :]

println("mean_Bii = ", mean(log.(fitted_all_filter.B0)),"\n", 
"var_Bii = ", var(log.(fitted_all_filter.B0))/abs(mean(log.(fitted_all_filter.B0))), "\n",
"mean_Eii = ", mean(fitted_all_filter.E), "\n",
"var_Eii = ", var(fitted_all_filter.E)/abs(mean(fitted_all_filter.E)), "\n",
"cor_ii = ", cor(log.(fitted_all_filter.B0), fitted_all_filter.E), "\n")

Nα, B_m, E_m, T_m, Ed,temp_all, allα = get_init_param(all_ℵii, num_temps)
E_p = mean(fitted_all_filter.E)
T_p = mean(fitted_all_filter.Th)
init_in = [exp(B_m), E_p, T_p, Ed]
fit_ij = curve_fit(temp_SS, temp_all, allα, init_in)
## calculate_r2
r_square = calculate_r2(fit_ij, temp_all, allα)
pred = abs.(temp_SS(temp, fit_ij.param))
E = fit_ij.param[2]

f = Figure(fontsize = 35,size = (1200, 1200));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|αij|)", ygridvisible = false, xgridvisible = false, xlabelsize = 50, ylabelsize = 50)
scatter!(ax, temp_all .- 273.15, log.(abs.(allα)), color = "#285C93", alpha = 0.5)
lines!(ax, Temp_rich, log.(pred), color = ("#E17542", 1), linewidth = 1)
text!(10, -4.5, text = "E = $(round(E, digits = 3))", align = (:center, :center), fontsize = 35)
f

###################
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

T = 273.15 + 25
ρ_t = [0.0000 0.0000]
Random.seed!(6)
p_0 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
prob = ODEProblem(dxx!, x0, tspan, p_0)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
sur_0 = (1:N)[bm .> 1.0e-7]
tempe_d_sur_0 = vcat(tempe_d[sur_0, sur_0]...)

ρ_t = [-0.9999 -0.9999]
Random.seed!(6)
p_1 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
prob = ODEProblem(dxx!, x0, tspan, p_1)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
sur_1 = (1:N)[bm .> 1.0e-7]
tempe_d_sur_1 = vcat(tempe_d[sur_1, sur_1]...)

path_ii = ["../results/αii_fitted_0.csv", "../results/αii_fitted-1.csv"]
path_ij = ["../results/αij_fitted_all_0.csv", "../results/αij_fitted_all-1.csv"]

fitted_0ii = filter_fitted(path_ii[1], p_0)
fitted_0ij = filter_fitted(path_ij[1], p_0)
fitted_0 = vcat(fitted_0ii[:,1:7], fitted_0ij)
fitted_0r = filter_fitted("../results/r_fitted0.csv", p_0)
fitted_1ii = filter_fitted(path_ii[2], p_1)
fitted_1ij = filter_fitted(path_ij[2], p_1)
fitted_1 = vcat(fitted_1ii[:,1:7], fitted_1ij)
fitted_1r = filter_fitted("../results/r_fitted-1.csv", p_1)

# E_u0 = vcat([fill(p_0.E[i,1], N) for i in 1:N]...)[fitted_0.id]
# E_u1 = vcat([fill(p_1.E[i,1], N) for i in 1:N]...)[fitted_1.id]
# # plot(E_u0, fitted_0.E )
# fitted_0[in.(fitted_0.id, Ref(tempe_d_sur_0)), :]
# fitted_1[in.(fitted_1.id, Ref(tempe_d_sur_1)), :].E

######### B0 ################
f = Figure(fontsize = 30, size = (1800, 600));
Label(f[:,0], "Minimal Trade-off", fontsize = 50, rotation = pi/2)
# Label(f[:,0, TopLeft()], "(a)")
ax1 = Axis(f[1,1], xlabel = "log(|B₀|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
xlims!(ax1, -20, nothing)
density!(ax1, log.(fitted_0ii.B0), label = L"α_{ii}", color = ("#FA8328", 0.8), strokewidth = 3, strokecolor = "#FA8328")
lines!(ax1, [log(median(fitted_0ii.B0)), log.(median(fitted_0ii.B0))], [0, 0.35], linestyle = :dash, color = ("#FA8328", 1.0), linewidth = 5)
lines!(ax1, [log(median(fitted_0ii.B0)), -3],[0.2, 0.2], linestyle = :dot, color = ("#FA8328", 0.9), linewidth = 3)
text!(ax1, -3, 0.2 , text = "$(round(log(median(fitted_0ii.B0)),digits = 2))", align = (:left, :center), fontsize = 20, color = "#FA8328")
density!(ax1, log.(fitted_0ij.B0), label = L"α_{i ≠ j}", color = ("#015845", 0.5), strokewidth = 3, strokecolor = "#015845")
lines!(ax1, [log.(median(fitted_0ij.B0)), log.(median(fitted_0ij.B0))], [0, 0.35], linestyle = :dash, color = ("#015845", 0.9), linewidth = 5)
lines!(ax1, [log(median(fitted_0ij.B0)), -3],[0.15, 0.15], linestyle = :dot, color = ("#015845", 0.9), linewidth = 3)
text!(ax1, -3, 0.15 , text = "$(round(log(median(fitted_0ij.B0)),digits = 2))", align = (:left, :center), fontsize = 20, color = "#015845")
axislegend(position = :rt)
ax2 = Axis(f[1,2], xlabel = "E", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax2, fitted_0ii.E, label = L"α_{ii}", color = ("#FA8328", 0.8), strokewidth = 3, strokecolor = "#FA8328")
lines!(ax2, [median(fitted_0ii.E), median(fitted_0ii.E)], [0, 0.7], linestyle = :dash, color = ("#FA8328", 1.0), linewidth = 5)
lines!(ax2, [median(fitted_0ii.E), 7],[0.5, 0.5], linestyle = :dot, color = ("#FA8328", 0.9), linewidth = 3)
text!(ax2, 7, 0.5 , text = "$(round(median(fitted_0ii.E),digits = 2)) eV", align = (:left, :center), fontsize = 20, color = "#FA8328")
density!(ax2, fitted_0ij.E, label = L"α_{i ≠ j}", color = ("#015845", 0.5), strokewidth = 3, strokecolor = "#015845")
lines!(ax2, [median(fitted_0ij.E), median(fitted_0ij.E)], [0, 0.7], linestyle = :dash, color = ("#015845", 0.9), linewidth = 5)
lines!(ax2, [median(fitted_0ij.E), 7],[0.4, 0.4], linestyle = :dot, color = ("#015845", 0.9), linewidth = 3)
text!(ax2, 7, 0.4, text = "$(round(median(fitted_0ij.E),digits = 2)) eV", align = (:left, :center), fontsize = 20, color = "#015845")
axislegend(position = :rt)
ax3 = Axis(f[1,3], xlabel = "Tₚₖ", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax3, fitted_0ii.Th .- 273.15, label = L"α_{ii}", color = ("#FA8328", 0.8), strokewidth = 3, strokecolor = "#FA8328")
lines!(ax3, [median(fitted_0ii.Th) - 273.15, median(fitted_0ii.Th) - 273.15], [0, 0.1], linestyle = :dash, color = ("#FA8328", 1.0), linewidth = 5)
lines!(ax3, [0, median(fitted_0ii.Th) - 273.15],[0.09, 0.09], linestyle = :dot, color = ("#FA8328", 0.9), linewidth = 3)
text!(ax3, 0, 0.09 , text = "$(round(median(fitted_0ii.Th) - 273.15,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#FA8328")
density!(ax3, fitted_0ij.Th .- 273.15, label = L"α_{i ≠ j}", color = ("#015845", 0.5), strokewidth = 3, strokecolor = "#015845")
lines!(ax3, [median(fitted_0ij.Th)- 273.15, median(fitted_0ij.Th)- 273.15], [0, 0.1], linestyle = :dash, color = ("#015845", 0.9), linewidth = 5)
lines!(ax3, [0, median(fitted_0ij.Th) - 273.15],[0.07, 0.07], linestyle = :dot, color = ("#015845", 0.9), linewidth = 3)
text!(ax3, 0, 0.07 , text = "$(round(median(fitted_0ij.Th) - 273.15,digits = 2)) °C", align = (:right, :center), fontsize = 20, color = "#015845")
axislegend(position = :rt)
f
save("../results/TPCα0_iiij.pdf", f) 


f = Figure(fontsize = 30, size = (1800, 600));
Label(f[:,0], "Maximal Trade-off", fontsize = 50, rotation = pi/2)
ax1 = Axis(f[1,1], xlabel = "log(|B₀|)", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax1, log.(fitted_1ii.B0), label = L"α_{ii}", color = ("#FA8328", 0.8), strokewidth = 3, strokecolor = "#FA8328")
lines!(ax1, [log(median(fitted_1ii.B0)), log.(median(fitted_1ii.B0))], [0, 0.8], linestyle = :dash, color = ("#FA8328", 1.0), linewidth = 5)
lines!(ax1, [-12, log(median(fitted_1ii.B0))],[0.6, 0.6], linestyle = :dot, color = ("#FA8328", 0.9), linewidth = 3)
text!(ax1, -12, 0.6 , text = "$(round(log(median(fitted_1ii.B0)),digits = 2))", align = (:right, :center), fontsize = 20, color = "#FA8328")
density!(ax1, log.(fitted_1ij.B0), label = L"α_{i ≠ j}", color = ("#015845", 0.5), strokewidth = 3, strokecolor = "#015845")
lines!(ax1, [log.(median(fitted_1ij.B0)), log.(median(fitted_1ij.B0))], [0, 0.8], linestyle = :dash, color = ("#015845", 0.9), linewidth = 5)
lines!(ax1, [-12, log(median(fitted_1ij.B0))],[0.5, 0.5], linestyle = :dot, color = ("#015845", 0.9), linewidth = 3)
text!(ax1, -12, 0.5 , text = "$(round(log(median(fitted_1ij.B0)),digits = 2))", align = (:right, :center), fontsize = 20, color = "#015845")
axislegend(position = :rt)
# Label(f[1,1, TopLeft()], "(a)")
ax2 = Axis(f[1,2], xlabel = "E", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax2, fitted_1ii.E, label = L"α_{ii}", color = ("#FA8328", 0.8), strokewidth = 3, strokecolor = "#FA8328")
lines!(ax2, [median(fitted_1ii.E), median(fitted_1ii.E)], [0, 0.8], linestyle = :dash, color = ("#FA8328", 1.0), linewidth = 5)
lines!(ax2, [median(fitted_1ii.E), 5],[0.5, 0.5], linestyle = :dot, color = ("#FA8328", 0.9), linewidth = 3)
text!(ax2, 5, 0.5 , text = "$(round(median(fitted_1ii.E),digits = 2)) eV", align = (:left, :center), fontsize = 20, color = "#FA8328")
density!(ax2, fitted_1ij.E, label = L"α_{i ≠ j}", color = ("#015845", 0.5), strokewidth = 3, strokecolor = "#015845")
lines!(ax2, [median(fitted_1ij.E), median(fitted_1ij.E)], [0, 0.8], linestyle = :dash, color = ("#015845", 0.9), linewidth = 5)
lines!(ax2, [median(fitted_1ij.E), 5],[0.4, 0.4], linestyle = :dot, color = ("#015845", 0.9), linewidth = 3)
text!(ax2, 5, 0.4, text = "$(round(median(fitted_1ij.E),digits = 2)) eV", align = (:left, :center), fontsize = 20, color = "#015845")
axislegend(position = :rt)
ax3 = Axis(f[1,3], xlabel = "Tₚₖ", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
density!(ax3, fitted_1ii.Th .- 273.15, label = L"α_{ii}", color = ("#FA8328", 0.8), strokewidth = 3, strokecolor = "#FA8328")
lines!(ax3, [median(fitted_1ii.Th)- 273.15, median(fitted_1ii.Th)- 273.15], [0, 0.17], linestyle = :dash, color = ("#FA8328", 1.0), linewidth = 5)
lines!(ax3, [median(fitted_1ii.Th) - 273.15, 40],[0.1, 0.1], linestyle = :dot, color = ("#FA8328", 0.9), linewidth = 3)
text!(ax3, 40, 0.1 , text = "$(round(median(fitted_1ii.Th) - 273.15,digits = 2)) °C", align = (:left, :center), fontsize = 20, color = "#FA8328")
density!(ax3, fitted_1ij.Th .- 273.15, label = L"α_{i ≠ j}", color = ("#015845", 0.5), strokewidth = 3, strokecolor = "#015845")
lines!(ax3, [median(fitted_1ij.Th)- 273.15, median(fitted_1ij.Th)- 273.15], [0, 0.17], linestyle = :dash, color = ("#015845", 0.9), linewidth = 5)
lines!(ax3, [median(fitted_1ij.Th) - 273.15, 40],[0.08, 0.08], linestyle = :dot, color = ("#015845", 0.9), linewidth = 3)
text!(ax3, 40, 0.08 , text = "$(round(median(fitted_1ij.Th) - 273.15,digits = 2)) °C", align = (:left, :center), fontsize = 20, color = "#015845")
axislegend(position = :rt)
f
save("../results/TPCα-1_iiij.pdf", f) 




##################

f = Figure(fontsize = 30, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "E", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
xlims!(nothing, 7)
# density!(ax1, fitted_0r.E, label = "Er", color = ("#20090E", 0.2), strokewidth = 3, strokecolor = "#20090E")
# lines!(ax1, [median(fitted_0r.E), median(fitted_0r.E)], [0, 1.6], linestyle = :dash, color = ("#20090E", 0.9), linewidth = 5)
density!(ax1, p_0.E[:,1], label = "Eu", color = ("#FF9776", 0.6), strokewidth = 3, strokecolor = "#FF9776")
lines!(ax1, [median(p_0.E[:,1]), median(p_0.E[:,1])], [0, 1.6], linestyle = :dash, color = ("#AD5525", 0.9), linewidth = 5)
density!(ax1, p_0.E[:,2], label = "Em", color = ("#E8C99E", 0.5), strokewidth = 3, strokecolor = "#E8C99E")
lines!(ax1, [median(p_0.E[:,2]), median(p_0.E[:,2])], [0, 1.6], linestyle = :dash, color = ("#F8BA17", 0.9), linewidth = 5)
density!(ax1, fitted_0.E, label = "Eα", color = ("#C25E8B", 0.8), strokewidth = 3, strokecolor = "#C25E8B")
lines!(ax1, [median(fitted_0.E), median(fitted_0.E)], [0, 1.6], linestyle = :dash, color = ("#601210", 1.0), linewidth = 5)
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/Eau0.pdf", f)
f = Figure(fontsize = 30, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "E", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
xlims!(nothing, 7)
density!(ax1, fitted_1.E, label = "Eα", color = ("#C25E8B", 0.3), strokewidth = 3, strokecolor = "#C25E8B")
lines!(ax1, [median(fitted_1.E), median(fitted_1.E)], [0, 1.5], linestyle = :dash, color = ("#601210", 1.0), linewidth = 5)
# density!(ax1, fitted_1r.E, label = "Er", color = ("#20090E", 0.2), strokewidth = 3, strokecolor = "#20090E")
# lines!(ax1, [median(fitted_1r.E), median(fitted_1r.E)], [0, 1.5], linestyle = :dash, color = ("#20090E", 0.9), linewidth = 5)
density!(ax1, p_1.E[:,1], label = "Eu", color = ("#FF9776", 0.6), strokewidth = 3, strokecolor = "#FF9776")
lines!(ax1, [median(p_1.E[:,1]), median(p_1.E[:,1])], [0, 1.5], linestyle = :dash, color = ("#AD5525", 0.9), linewidth = 5)
density!(ax1, p_1.E[:,2], label = "Em", color = ("#E8C99E", 0.5), strokewidth = 3, strokecolor = "#E8C99E")
lines!(ax1, [median(p_1.E[:,2]), median(p_1.E[:,2])], [0, 1.5], linestyle = :dash, color = ("#F8BA17", 0.9), linewidth = 5)
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/Eau-1.pdf", f)

f = Figure(fontsize = 30, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Tₚₖ", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
# xlims!(nothing, 7)
density!(ax1, p_0.Tp[:,1] .- 273.15, label = "Tpu", color = ("#82AC6D", 0.5), strokewidth = 3, strokecolor = "#82AC6D")
lines!(ax1, [median(p_0.Tp[:,1])- 273.15, median(p_0.Tp[:,1])- 273.15], [0, 0.1], linestyle = :dash, color = ("#12473D", 0.9), linewidth = 5)
density!(ax1, p_0.Tp[:,2] .- 273.15, label = "Tpm", color = ("#C1C6E8", 0.6), strokewidth = 3, strokecolor = "#C1C6E8")
lines!(ax1, [median(p_0.Tp[:,2])- 273.15, median(p_0.Tp[:,2])- 273.15], [0, 0.1], linestyle = :dash, color = ("#9585B4", 0.9), linewidth = 5)
density!(ax1, fitted_0.Th .- 273.15, label = "Tpα", color = ("#5676A5", 0.8), strokewidth = 3, strokecolor = "#5676A5")
lines!(ax1, [median(fitted_0.Th)- 273.15, median(fitted_0.Th)- 273.15], [0, 0.1], linestyle = :dash, color = ("#0758AE", 0.9), linewidth = 5)
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(c)")
f
save("../results/Tpau0.pdf", f)
f = Figure(fontsize = 30, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Tₚₖ", ylabel = "Density", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
# xlims!(nothing, 7)
density!(ax1, fitted_1.Th .- 273.15, label = "Tpα", color = ("#5676A5", 0.5), strokewidth = 3, strokecolor = "#5676A5")
lines!(ax1, [median(fitted_1.Th)- 273.15, median(fitted_1.Th)- 273.15], [0, 0.17], linestyle = :dash, color = ("#0758AE", 1.0), linewidth = 5)
density!(ax1, p_1.Tp[:,1] .- 273.15, label = "Tpu", color = ("#82AC6D", 0.5), strokewidth = 3, strokecolor = "#82AC6D")
lines!(ax1, [median(p_1.Tp[:,1])- 273.15, median(p_1.Tp[:,1])- 273.15], [0, 0.17], linestyle = :dash, color = ("#12473D", 0.9), linewidth = 5)
density!(ax1, p_1.Tp[:,2] .- 273.15, label = "Tpm", color = ("#C1C6E8", 0.6), strokewidth = 3, strokecolor = "#C1C6E8")
lines!(ax1, [median(p_1.Tp[:,2])- 273.15, median(p_1.Tp[:,2])- 273.15], [0, 0.17], linestyle = :dash, color = ("#9585B4", 0.9), linewidth = 5)
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(d)")
f
save("../results/Tpau-1.pdf", f)



f = Figure(fontsize = 30, size = (1800, 900));
ax1 = Axis(f[1,1], xlabel = "log(B₀)", ylabel = "frequency", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,2], xlabel = "E", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax3 = Axis(f[1,3], xlabel = "Tₚₖ", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax4 = Axis(f[1,1], xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
ax5 = Axis(f[1,2], xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
ax6 = Axis(f[1,3], xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)

hist!(ax1, log.(fitted_0ij.B0), color = ("#82AC6D", 0.7), bins = 100, label = "αᵢⱼ")
hist!(ax4, log.(fitted_0ii.B0), color = ("#C25E8B", 0.7), bins = 100, label = "αᵢᵢ")
hist!(ax2, fitted_0ij.E, color = ("#82AC6D", 0.7), bins = 100, label = "αᵢⱼ")
hist!(ax5, fitted_0ii.E, color = ("#C25E8B", 0.7), bins = 100, label = "αᵢᵢ")
hist!(ax3, fitted_0ij.Th .- 273.15, color = ("#82AC6D", 0.7), bins = 100, label = "αᵢⱼ")
hist!(ax6, fitted_0ii.Th .- 273.15, color = ("#C25E8B", 0.7), bins = 100, label = "αᵢᵢ")
axislegend(position = :rt)
f
save("../results/TPC_αij.pdf", f) 

hist([p.E[i,1] + p.E[j,1] for j in 1:N for i in 1:N])

fitted_0.E
all_E = [p.E[i,1] + p.E[j,1] for j in 1:N for i in 1:N]

f = Figure(fontsize = 30, size = (1800, 900));
ax1 = Axis(f[1,1], xlabel = "log(B0)", ylabel = "frequency", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,2], xlabel = "E", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax3 = Axis(f[1,3], xlabel = "Th", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
hist!(ax1, log.(fitted_0ii.B0), color = ("#82AC6D", 0.7), bins = 100, label = "ρ=0")
hist!(ax1, log.(fitted_1ii.B0), color = ("#C25E8B", 0.7), bins = 100, label = "ρ=-1")
hist!(ax2, fitted_0ii.E, color = ("#82AC6D", 0.7), bins = 100, label = "ρ=0")
hist!(ax2, fitted_1ii.E, color = ("#C25E8B", 0.7), bins = 100, label = "ρ=-1")
hist!(ax3, fitted_0ii.Th .- 273.15, color = ("#82AC6D", 0.7), bins = 100, label = "ρ=0")
hist!(ax3, fitted_1ii.Th .- 273.15, color = ("#C25E8B", 0.7), bins = 100, label = "ρ=-1")
axislegend(position = :rt)
f
save("../results/TPC_αii.pdf", f) 

################################
@load "../data/1com_0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;

num_temps = 38
Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)

allαii = abs.(vcat(all_ℵii...))
allαij = abs.(vcat(all_ℵij...))
temp_ii = vcat([repeat([temp[t]], length(all_ℵii[t])) for t in 1:num_temps]...)
temp_ij = vcat([repeat([temp[t]], length(all_ℵij[t])) for t in 1:num_temps]...)

f = Figure(fontsize = 35,size = (1200, 1200));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|α|)", ygridvisible = false, xgridvisible = false, xlabelsize = 50, ylabelsize = 50)
scatter!(ax, temp_ij .- 273.15, log.(abs.(allαij)), color = ("#015845", 0.4), label = "αij")
scatter!(ax, temp_ii .- 273.15, log.(abs.(allαii)), color = ("#FA8328", 0.5), label = "αii")
s1 = [MarkerElement(color = ("#FA8328", 0.8), marker = :circle)]
s2 = [MarkerElement(color = ("#015845", 0.8), marker = :circle)]
Legend(f[1,1], [s1, s2], tellheight = false, tellwidth = false, [ "αii", "αij"], halign = :left, valign = :top)
f
save("../results/TPCα0.pdf", f) 


fitted_all = CSV.read("../results/αij_fitted_all-1.csv", DataFrame, header=false)
df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_all = DataFrame(fitted_all, df_names);

fitted_all_filter = fitted_all[fitted_all.B0 .> 0, :]
fitted_all_filter = fitted_all_filter[fitted_all_filter.E .> 0, :]
fitted_all_filter = fitted_all_filter[fitted_all_filter.Ed .> 0, :]
fitted_all_filter = fitted_all_filter[fitted_all_filter.r2 .> 0.9, :]

fitted_all_filter[fitted_all_filter.Th .> 0, :]
