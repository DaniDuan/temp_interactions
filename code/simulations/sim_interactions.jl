include("./sim_frame.jl");
using ProgressMeter, RCall

N=100
M=50
L = 0.3
### Temp params 
num_temps = 38
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)


all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); 
all_up_ℵij = Vector{Vector{Float64}}(); all_low_ℵij = Vector{Vector{Float64}}(); 
all_ℵij_sum = Vector{Vector{Float64}}(); all_D_ℵij = Vector{Vector{Float64}}();
all_ℵii_sur = Vector{Vector{Float64}}(); all_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); 
all_up_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); all_low_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); 
all_ℵij_sum_sur = Vector{Vector{Union{Float64, Missing}}}(); all_D_ℵij_sur =  Vector{Vector{Union{Float64, Missing}}}()
progress = Progress(num_temps; desc="Progress running:")
for i in range(0, stop = num_temps-1, length = num_temps)
    T = 273.15 + i
    Random.seed!(6)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    # println(bm .> 1.0e-7)
    sur = (1:N)[bm .> 1.0e-7]
    N_s = length(sur)
    p_lv = Eff_LV_params(p=p, sol=sol);
    # number of species with r>0 at equilibium 
    ℵ = p_lv.ℵ
    ℵii = diag(ℵ)
    ℵij_sum = vec(sum(ℵ, dims = 2).- diag(ℵ))
    ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i != j]
    up_ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i > j]
    low_ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i < j]
    D_ℵij = up_ℵij./low_ℵij

    sur_ℵ = ℵ[sur, sur]
    ℵii_sur = diag(sur_ℵ)
    if N_s > 1
        ℵij_sum_sur = vec(sum(sur_ℵ, dims = 2).- diag(sur_ℵ))
        ℵij_sur = [sur_ℵ[i,j] for i in 1:N_s for j in 1:N_s if i != j]
        up_ℵij_sur = [sur_ℵ[i,j] for i in 1:N_s for j in 1:N_s if i > j]
        low_ℵij_sur = [sur_ℵ[i,j] for i in 1:N_s for j in 1:N_s if i < j]
        D_ℵij_sur = up_ℵij_sur./low_ℵij_sur
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_up_ℵij, up_ℵij); push!(all_low_ℵij, low_ℵij); 
        push!(all_ℵij_sum, ℵij_sum); push!(all_D_ℵij, D_ℵij);
        push!(all_ℵii_sur, ℵii_sur); push!(all_ℵij_sur, ℵij_sur); push!(all_up_ℵij_sur, up_ℵij_sur); push!(all_low_ℵij_sur, low_ℵij_sur); 
        push!(all_ℵij_sum_sur, ℵij_sum_sur); push!(all_D_ℵij_sur, D_ℵij_sur);
    else 
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_up_ℵij, up_ℵij); push!(all_low_ℵij, low_ℵij); 
        push!(all_ℵij_sum, ℵij_sum); push!(all_D_ℵij, D_ℵij);
        push!(all_ℵii_sur, ℵii_sur); push!(all_ℵij_sur, [missing]); push!(all_up_ℵij_sur, [missing]); push!(all_low_ℵij_sur, [missing]); 
        push!(all_ℵij_sum_sur, [missing]); push!(all_D_ℵij_sur, [missing]);
    end
    next!(progress)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"


@load "../data/1com_0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur
D = (all_ℵii = all_ℵii, all_ℵij = all_ℵij, all_up_ℵij = all_up_ℵij,all_low_ℵij = all_low_ℵij, all_ℵij_sum= all_ℵij_sum, all_D_ℵij = all_D_ℵij,
    all_ℵii_sur = all_ℵii_sur,  all_ℵij_sur = all_ℵij_sur, all_up_ℵij_sur = all_up_ℵij_sur, all_low_ℵij_sur = all_low_ℵij_sur, all_ℵij_sum_sur = all_ℵij_sum_sur, all_D_ℵij_sur = all_D_ℵij_sur);
Dnames = ("αii", "αij", "up_αij", "low_αij", "sum_αij", "up_low", "αii_sur", "αij_sur", "up_αij_sur", "low_αij_sur", "sum_αij_sur", "up_low_sur");

### @save "../data/1com-1.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur

Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)

########### Fitting community α #################
include("./fitting.jl");

all_ℵij_sum
# f1
progress = Progress(N; desc="Progress running:")
temp = collect(Temp_rich .+273.15)
f1 = Figure(size = (1200, 1200));
fitted = zeros(Float64, N, 6)
@time for i in 1:N 
        next!(progress)
        αii = [all_ℵij_sum[t][i] for t in 1:num_temps]
        Nα, init_in, AIC_in, temp_all, allα = try_params(αii, num_temps, 2000)
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
        fitted[Int(i),:] = vcat(params, aic_value, r_square)
        ## plotting data
        ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
        scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
        lines!(ax1, Temp_rich, pred, color = ("#E17542", 1), linewidth = 1)
    end 
R"library(beepr); beep(sound = 4, expr = NULL)"

df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted = DataFrame(fitted, df_names);
# CSV.write("../results/αii_fitted_0.csv", fitted, writeheader=false)
CSV.write("../results/αij_sum_fitted_0.csv", fitted, writeheader=false)

f1
# save("/Users/Danica/Documents/temp_interactions/results/αii_fitted_0.png", f1) 
save("/Users/Danica/Documents/temp_interactions/results/αij_sum_fitted_0.png", f1) 
# fitted = CSV.read("../results/αii_fitted-1.csv", DataFrame, header=false)

fitted = fitted[fitted.E .> 0, :]
fitted = fitted[fitted.B0 .> 0, :]

mean_Bii = mean(log.(fitted.B0))
var_Bii = var(log.(fitted.B0))/abs(mean(log.(fitted.B0)))
mean_Eii = mean(fitted.E)
var_Eii = var(fitted.E)/abs(mean(fitted.E))
cor_ii = cor(log.(fitted.B0), fitted.E)

Nα, B_m, E_m, T_m, Ed,temp_all, allα = get_init_param(all_ℵii, num_temps)
E_p = mean(fitted.E[fitted.E .> 0])
T_p = mean(fitted.Th[fitted.Th .> 0])
init_in = [exp(B_m), E_p, T_p, Ed]
fit_ii = curve_fit(temp_SS, temp_all, allα, init_in)
## calculate_r2
r_square = calculate_r2(fit_ii, temp_all, allα)
pred = abs.(temp_SS(temp, fit_ii.param))
E = fit_ii.param[2]

f = Figure(fontsize = 35,size = (1200, 1200));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|αii|)", ygridvisible = false, xgridvisible = false, xlabelsize = 50, ylabelsize = 50)
scatter!(ax, temp_all .- 273.15, log.(abs.(allα)), color = "#285C93", alpha = 0.5)
lines!(ax, Temp_rich, log.(pred), color = ("#E17542", 1), linewidth = 1)
text!(10, -4.5, text = "E = $(round(E, digits = 3))", align = (:center, :center), fontsize = 35)
f
# save("/Users/Danica/Documents/temp_interactions/results/αii_fitted_all-1.png", f) 
