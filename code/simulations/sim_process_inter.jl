include("./fitting.jl");
using ProgressMeter, RCall

N=100
M=50
L = 0.3
### Temp params 
num_temps = 38
Tr=273.15+10; Ed=3.5 

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
function filter_fitted(path, p)
    fitted = CSV.read(path, DataFrame, header=false)
    df_names = ["B0","E","Th","Ed","AIC","r2"]
    fitted = rename(fitted, df_names);
    if nrow(fitted) < 101
        fitted[!,:"Bu"] = p.B[:,1]
        fitted[!,:"Eu"] = p.E[:,1]
        fitted[!,:"Thu"] = p.Tp[:,1]
    end
    fitted = fitted[fitted.E .> eps(), :]
    fitted = fitted[fitted.B0 .> eps(), :]
    fitted = fitted[fitted.Th .< 323.15, :]
    fitted = fitted[fitted.Th .> 273.15, :]
    fitted = fitted[fitted.r2 .> 0.90, :]
    return fitted
end

ρ_t = [0.0000 0.0000]
Random.seed!(6)
p_0 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
ρ_t = [-0.9999 -0.9999]
Random.seed!(6)
p_1 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)

path_ii = ["../results/αii_fitted_0.csv", "../results/αii_fitted-1.csv"]
path_ij = ["../results/αij_fitted_all_0_new.csv", "../results/αij_fitted_all-1.csv"]

fitted_0ii = filter_fitted(path_ii[1], p_0)
fitted_0ij = filter_fitted(path_ij[1], p_0)
fitted_0 = vcat(fitted_0ii[:,1:6], fitted_0ij)
fitted_1ii = filter_fitted(path_ii[2], p_1)
fitted_1ij = filter_fitted(path_ij[2], p_1)
fitted_1 = vcat(fitted_1ii[:,1:6], fitted_1ij)

# plot(fitted_0ii.E, fitted_0ii.Eu)

f = Figure(fontsize = 30, size = (1800, 900));
ax1 = Axis(f[1,1], xlabel = "log(B0)", ylabel = "frequency", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,2], xlabel = "E", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax3 = Axis(f[1,3], xlabel = "Th", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
hist!(ax1, log.(fitted_0.B0), color = ("#82AC6D", 0.7), bins = 100, label = "ρ=0")
hist!(ax1, log.(fitted_1.B0), color = ("#C25E8B", 0.7), bins = 100, label = "ρ=-1")
hist!(ax2, fitted_0.E, color = ("#82AC6D", 0.7), bins = 100, label = "ρ=0")
hist!(ax2, fitted_1.E, color = ("#C25E8B", 0.7), bins = 100, label = "ρ=-1")
hist!(ax3, fitted_0.Th .- 273.15, color = ("#82AC6D", 0.7), bins = 100, label = "ρ=0")
hist!(ax3, fitted_1.Th .- 273.15, color = ("#C25E8B", 0.7), bins = 100, label = "ρ=-1")
axislegend(position = :rt)
f
save("../results/TPC_α.pdf", f) 

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


########################################
num_temps = 38
Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)

@load "../data/1com_0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
D = (all_ℵii = all_ℵii, all_ℵij = all_ℵij, all_up_ℵij = all_up_ℵij,all_low_ℵij = all_low_ℵij, all_ℵij_sum= all_ℵij_sum, all_D_ℵij = all_D_ℵij,
    all_ℵii_sur = all_ℵii_sur,  all_ℵij_sur = all_ℵij_sur, all_up_ℵij_sur = all_up_ℵij_sur, all_low_ℵij_sur = all_low_ℵij_sur, all_ℵij_sum_sur = all_ℵij_sum_sur, all_D_ℵij_sur = all_D_ℵij_sur);
Dnames = ("αii", "αij", "up_αij", "low_αij", "sum_αij", "up_low", "αii_sur", "αij_sur", "up_αij_sur", "low_αij_sur", "sum_αij_sur", "up_low_sur");

all_ℵ = [vcat(all_ℵii[t], all_ℵij[t]) for t in 1: num_temps]
fitted_ii = CSV.read("../results/αii_fitted_0.csv", DataFrame, header=false)
fitted_ij = CSV.read("../results/αij_fitted_all_0.csv", DataFrame, header=false)
fitted = vcat(fitted_ii, fitted_ij)
df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted = DataFrame(fitted, df_names);
temp = collect(Temp_rich .+273.15)
f1 = Figure(size = (800, 800));
for i in 1:25
    # Random.seed!(i*5)
    Random.seed!(i-1)
    n = rand(1:10000)
    params = fitted[Int(n),1:4]
    α = [all_ℵ[t][n] for t in 1:num_temps]
    pred = abs.(temp_SS(temp, params))
    ax1 = Axis(f1[Int(floor((i-1)/5+1)),Int((i-1) % 5+1)], ygridvisible = false, xgridvisible = false)
    scatter!(ax1, Temp_rich, abs.(α), color = "#5676A5", alpha = 0.5)
    lines!(ax1, Temp_rich, pred, color = ("#C25E8B", 1), linewidth = 2)
end 
Label(f1[1,1, TopLeft()], "(a)", fontsize = 25)

f1

save("../results/TPCα0_example.pdf", f1) 

