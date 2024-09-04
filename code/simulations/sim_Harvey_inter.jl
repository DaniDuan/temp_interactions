include("./fitting.jl");
# using ProgressMeter
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
###################################
# Retrieve the environment variable as a string
index_str = ENV["SLURM_ARRAY_TASK_ID"]
# Convert the string to a numeric value (e.g., Integer)
index = parse(Int, index_str)

@load "../data/1com_0.jld2" all_ℵij
# f = Figure()
# ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "upper/lower")
# Nα, B_m, E_m, T_m, Ed,temp_all, allα = get_init_param(all_ℵij, num_temps)
# scatter!(ax, temp_all.-273.15, log.(abs.(allα)), color = "#285C93", alpha = 0.5)
# f

# progress = Progress(10; desc="Progress running:")
# f1 = Figure(size = (1200, 1200));

fitted = zeros(Float64, 10, 6)
# for index in 1: 990
for n in 1:10
    # next!(progress)
    i = 10*(index-1) + n
    αii = [all_ℵij[t][i] for t in 1:num_temps]
    Nα, init_in, AIC_in, temp_all, allα = try_params(αii, num_temps, 2000)
    try
        fit_ii = curve_fit(temp_SS, temp_all, allα, init_in)
        r_square = calculate_r2(fit_ii, temp_all, allα)
        params = fit_ii.param
    ## calculate_r2
    pred = abs.(temp_SS(temp, params))
    ss_res = sum((allα .- pred).^2)
    ss_tot = sum((allα .- mean(allα)).^2)
    r_square = 1 - ss_res / ss_tot
    ## calculate_AIC
    aic_value = N * log(ss_res / N) + 2 * 4
    ## store data
    fitted[Int(n),:] = vcat(params, aic_value, r_square)
    catch e
        fitted[Int(n),:] = [0.0, 0.0, 0.0, 0.0, Inf, 0.0]
    end 
    ## plotting data
    # ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    # scatter!(ax1, temp_all, abs.(allα), color = "#285C93", alpha = 0.5)
    # lines!(ax1, Temp_rich, pred, color = ("#E17542", 1), linewidth = 1)
end 

@save "../results/20240902/inter_ij_0/fit_allij_0_$(index).jld2" fitted
# df_names = ["B0","E","Th","Ed","AIC","r2"]
# fitted = DataFrame(fitted, df_names);
# CSV.write("../results/αij_fitted-1.csv", fitted, writeheader=false)
