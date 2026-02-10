include("./fitting.jl");
using ProgressMeter, LsqFit, GLM
N=100;
M=50;
L = 0.3;
### Temp params 
num_temps = 38;
Tr=273.15+10; Ed=3.5;
Temp_rich = range(0, num_temps-1, length = num_temps);
k = 0.0000862 # Boltzman constant
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr);
temp = collect(Temp_rich .+273.15);
temp_SS(T, params) = params[1] .* exp.((-params[2]./k) * ((1 ./T) .-(1/Tr)))./(1 .+ (params[2]./(params[4] .- params[2])) .* exp.(params[4]/k * (1 ./params[3] .- 1 ./T)));

@load "../data/1com_0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
D = (all_ℵii = all_ℵii, all_ℵij = all_ℵij, all_up_ℵij = all_up_ℵij,all_low_ℵij = all_low_ℵij, all_ℵij_sum= all_ℵij_sum, all_D_ℵij = all_D_ℵij,
    all_ℵii_sur = all_ℵii_sur,  all_ℵij_sur = all_ℵij_sur, all_up_ℵij_sur = all_up_ℵij_sur, all_low_ℵij_sur = all_low_ℵij_sur, all_ℵij_sum_sur = all_ℵij_sum_sur, all_D_ℵij_sur = all_D_ℵij_sur);
Dnames = ("αii", "αij", "up_αij", "low_αij", "sum_αij", "up_low", "αii_sur", "αij_sur", "up_αij_sur", "low_αij_sur", "sum_αij_sur", "up_low_sur");

f1 = Figure(size = (1200, 1200));
fitted = CSV.read("../results/αii_fitted_0.csv", DataFrame, header=false)
temp = collect(Temp_rich .+273.15)
for i in 1:N  
    params = fitted[Int(i),1:4]
    αii = [all_ℵii[t][i] for t in 1:num_temps]
    pred = abs.(temp_SS(temp, params))
    ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
    lines!(ax1, Temp_rich, pred, color = ("#E17542", 1), linewidth = 2)
end 
display(f1);
save("../results/TPCαii0.pdf", f1) 

f2 = Figure(size = (1200, 1200));
fitted = CSV.read("../results/αij_sum_fitted_0.csv", DataFrame, header=false)
temp = collect(Temp_rich .+273.15)
for i in 1:N 
    params = fitted[Int(i),1:4]
    αij = [all_ℵij_sum[t][i] for t in 1:num_temps]
    pred = abs.(temp_SS(temp, params))
    ax2 = Axis(f2[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    scatter!(ax2, Temp_rich, abs.(αij), color = "#285C93", alpha = 0.5)
    lines!(ax2, Temp_rich, pred, color = ("#015845", 1), linewidth = 2)
end 
display(f2);
save("../results/TPCαij0_sum.pdf", f2) 


############################################################
@load "../data/1com-1.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur;
D = (all_ℵii = all_ℵii, all_ℵij = all_ℵij, all_up_ℵij = all_up_ℵij,all_low_ℵij = all_low_ℵij, all_ℵij_sum= all_ℵij_sum, all_D_ℵij = all_D_ℵij,
    all_ℵii_sur = all_ℵii_sur,  all_ℵij_sur = all_ℵij_sur, all_up_ℵij_sur = all_up_ℵij_sur, all_low_ℵij_sur = all_low_ℵij_sur, all_ℵij_sum_sur = all_ℵij_sum_sur, all_D_ℵij_sur = all_D_ℵij_sur);
Dnames = ("αii", "αij", "up_αij", "low_αij", "sum_αij", "up_low", "αii_sur", "αij_sur", "up_αij_sur", "low_αij_sur", "sum_αij_sur", "up_low_sur");

f1 = Figure(size = (1200, 1200));
fitted = CSV.read("../results/αii_fitted-1.csv", DataFrame, header=false)
temp = collect(Temp_rich .+273.15)
for i in 1:N 
    params = fitted[Int(i),1:4]
    αii = [all_ℵii[t][i] for t in 1:num_temps]
    pred = abs.(temp_SS(temp, params))
    ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
    lines!(ax1, Temp_rich, pred, color = ("#E17542", 1), linewidth = 2)
end 
display(f1);
save("../results/TPCαii-1.pdf", f1) 

f2 = Figure(size = (1200, 1200));
fitted = CSV.read("../results/αij_sum_fitted-1.csv", DataFrame, header=false)
temp = collect(Temp_rich .+273.15)
for i in 1:N 
    params = fitted[Int(i),1:4]
    αij = [all_ℵij_sum[t][i] for t in 1:num_temps]
    pred = abs.(temp_SS(temp, params))
    ax2 = Axis(f2[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    scatter!(ax2, Temp_rich, abs.(αij), color = "#285C93", alpha = 0.5)
    lines!(ax2, Temp_rich, pred, color = ("#015845", 1), linewidth = 2)
end 
display(f2);
save("../results/TPCαij_sum-1.pdf", f2) 
