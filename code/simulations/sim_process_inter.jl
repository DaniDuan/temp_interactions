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
    @load "../results/fit_allij_0/fit_allij0_$(index).jld2" fitted
    i = 10*(index-1) + 10
    fitted_all[10*(index-1)+1:i,:] = fitted
end 

df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted_all = DataFrame(fitted_all, df_names);
CSV.write("../results/αij_fitted_all_0.csv", fitted_all, writeheader=false)

# fitted_all = CSV.read("../results/αij_fitted_all_re.csv", DataFrame, header=false)
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
