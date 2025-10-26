include("./fitting.jl");
using ProgressMeter, RCall
using ColorSchemes

N=100
M=50
# L = fill(0.3, N)
L = fill(0.0, N)
### Temp params 
num_temps = 38
Temp_rich = range(0, num_temps-1, length = num_temps)
# ρ_t= [0.0000 0.0000] #[-0.3500 -0.3500]; # realistic covariance
# ρ_t= [-0.9999 -0.9999]
ρ_t= [-0.5000 -0.5000]
niche = fill(1.0, M, N)
input_type = ["constant", "leaching", "chemostat", "self-renewing"]

Tr=273.15+10; Ed=3.5 
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); 
progress = Progress(num_temps; desc="Progress running:")
for i in range(0, stop = num_temps-1, length = num_temps)
    T = 273.15 + i
    Random.seed!(6)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, niche = niche, input_type = input_type[4])
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    sur = (1:N)[bm .> 1.0e-7]

    # println(bm .> 1.0e-7)
    p_lv = Eff_LV_params(p=p, sol=sol);
    # number of species with r>0 at equilibium 
    ℵ = p_lv.ℵ
    ℵii = diag(ℵ)
    ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i != j]

    push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); 
    next!(progress)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

@save "../data/1com_0_nol_rho4.jld2" all_ℵii all_ℵij

p_0 = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=273.15 + 25, ρ_t=[0.0000 0.0000], Tr=Tr, Ed=Ed)

progress = Progress(N; desc="Progress running:")
temp = collect(Temp_rich .+273.15)
f1 = Figure(size = (1200, 1200));
fitted = zeros(Float64, N, 6)
@time for i in 1:N 
    next!(progress)
    αii = [all_ℵii[t][i] for t in 1:num_temps]
    Nα, init_in, AIC_in, temp_all, allα = try_params(αii, num_temps, 2000)
    try
        fit_ii = curve_fit(temp_SS, temp_all, allα, init_in)
        r_square = calculate_r2(fit_ii, temp_all, allα)  # Correct R-squared calculation
        params = fit_ii.param
        
        ## calculate_AIC using correct residuals and sample size
        pred_at_data = temp_SS(temp_all, params)  # Predictions at actual data points
        ss_res = sum((allα .- pred_at_data).^2)
        aic_value = Nα * log(ss_res / Nα) + 2 * 4  # Use Nα instead of N
        
        ## store data
        fitted[Int(i),:] = vcat(params, aic_value, r_square)
        
        pred_full = abs.(temp_SS(temp, params))  # For plotting only
        ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
        scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
        lines!(ax1, Temp_rich, pred_full, color = ("#E17542", 1), linewidth = 1)
        
    catch e
        fitted[Int(i),:] = [0.0, 0.0, 0.0, 0.0, Inf, 0.0]
        
        ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
        scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
    end
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted = DataFrame(fitted, df_names);
CSV.write("../results/αii_fitted_nol_rho4.csv", fitted, writeheader=false)

f1
save("../results/αii_fitted_nol_rho4.pdf", f1) 

# @load "../data/1com_0_nol_rho4.jld2" all_ℵii all_ℵij
# all_ℵii_nol = all_ℵii; all_ℵij_nol = all_ℵij
# mean_ii_nol = [mean(all_ℵii_nol[i]) for i in 1:num_temps]
# mean_ij_nol = [mean(all_ℵij_nol[i]) for i in 1:num_temps]

# idx = collect(CartesianIndices(zeros(Float64, N, N)))
# ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
# dℵij_collect_nol = Vector{Vector{Float64}}() ; # all_R_collect_nol = Vector{Vector{Float64}}();
# @time for j in 1: num_temps
#     dℵij_H = Float64[]; #all_R_H = Float64[];
#     for i in 1:num_temps
#         all_ℵii = all_ℵii_nol[j]
#         all_ℵij = all_ℵij_nol[j]
#         A = zeros(Float64, N, N)
#         A[ind_off] = all_ℵij
#         A[diagind(A)] = all_ℵii
#         dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
#         append!(dℵij_H, dℵij);# append!(all_R_H, all_R[j]); 
#         next!(progress)
#     end 
#     push!(dℵij_collect_nol, dℵij_H);# push!(all_R_collect_0, all_R_H);
# end 
# R"library(beepr); beep(sound = 4, expr = NULL)"
# mean_dℵij_nol = [mean(dℵij_collect_nol[i]) for i in 1:num_temps]

# # @load "../data/1com_0.jld2" all_ℵii all_ℵij
# all_ℵii_l = all_ℵii; all_ℵij_l = all_ℵij
# mean_ii_l = [mean(all_ℵii_l[i]) for i in 1:num_temps]
# mean_ij_l = [mean(all_ℵij_l[i]) for i in 1:num_temps]
# dℵij_collect_l = Vector{Vector{Float64}}() ; # all_R_collect_l = Vector{Vector{Float64}}();
# @time for j in 1: num_temps
#     dℵij_H = Float64[]; #all_R_H = Float64[];
#     for i in 1:num_temps
#         all_ℵii = all_ℵii_l[j]
#         all_ℵij = all_ℵij_l[j]
#         A = zeros(Float64, N, N)
#         A[ind_off] = all_ℵij
#         A[diagind(A)] = all_ℵii
#         dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
#         append!(dℵij_H, dℵij);# append!(all_R_H, all_R[j]); 
#         next!(progress)
#     end 
#     push!(dℵij_collect_l, dℵij_H);# push!(all_R_collect_0, all_R_H);
# end 
# R"library(beepr); beep(sound = 4, expr = NULL)"
# mean_dℵij_l = [mean(dℵij_collect_l[i]) for i in 1:num_temps]

# f = Figure(fontsize = 35, resolution = (1200, 900));
# ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
# lines!(ax1, Temp_rich, mean_ij_nol, color = ("#285C93",1), linewidth = 5, label = "nol")
# lines!(ax1, Temp_rich, mean_ij_l, color = ("#E17542",1), linewidth = 5, label = "l")
# # axislegend(position = :rb)
# f
