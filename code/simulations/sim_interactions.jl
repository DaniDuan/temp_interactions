include("./sim_frame.jl");
using ProgressMeter, RCall

N=100
M=50
L = 0.3
### Temp params 
num_temps = 38
ρ_t= [-0.9999 -0.9999]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
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


# @load "../data/1com-1.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur
using JSON
D = (all_ℵii = all_ℵii, all_ℵij = all_ℵij, all_up_ℵij = all_up_ℵij,all_low_ℵij = all_low_ℵij, all_ℵij_sum= all_ℵij_sum, all_D_ℵij = all_D_ℵij,
    all_ℵii_sur = all_ℵii_sur,  all_ℵij_sur = all_ℵij_sur, all_up_ℵij_sur = all_up_ℵij_sur, all_low_ℵij_sur = all_low_ℵij_sur, all_ℵij_sum_sur = all_ℵij_sum_sur, all_D_ℵij_sur = all_D_ℵij_sur);
Dnames = ("αii", "αij", "up_αij", "low_αij", "sum_αij", "up_low", "αii_sur", "αij_sur", "up_αij_sur", "low_αij_sur", "sum_αij_sur", "up_low_sur");

### saving data in different formats 
# using JSON
# json_d = JSON.json(D)
# open("../data/1com-1.json", "w") do file 
#     write(file, json_d)
# end 
# # @save "../data/1com-1.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur

Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)

########### Fitting community α #################
include("./fitting.jl");

# iters_v = 500
# progress = Progress(length(D); desc="running progress:")
# for fn in 1:length(D)
#     Nα, init_in, AIC_in, temp_all, allα = try_params(D[fn], num_temps, iters_v)
#     name = Symbol("param_", Dnames[fn])
#     fit_ss_final = curve_fit(temp_SS, temp_all, allα, init_in)
#     name_final = Symbol("result_", Dnames[fn])
#     @eval $name_final = fit_ss_final.param
#     next!(progress)
# end 

# f1 = Figure(resolution = (1200, 1200));
# ax1 = Axis(f1[1,1], xlabel = "Temperature (°C)", ylabel = "αii", xlabelsize = 50, ylabelsize = 50)
# scatter!(ax1, temp_all .- 273.15, log.(abs.(allii)), color = "#285C93", alpha = 0.5)
# lines!(ax1, Temp_rich, log.(abs.(temp_SS(temp, fit_ii.param))), color = ("#E17542", 1), linewidth = 1)
# lines!(ax1, Temp_rich, ii, color = ("#285C93", 1), linewidth = 1)
# f1

################# Fitting each α ##################
# progress = Progress(length(D); desc="Progress running:")
# for fn in 1:length(D)
#     f = Figure(size = (1200, 1200));
#     ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "$(Dnames[fn])", xlabelsize = 50, ylabelsize = 50)
#     Nα, B_m, E_up, T_m, temp_all, allα = get_init_param(D[fn], num_temps)
#         # ax = Axis(f[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
#     scatter!(ax, temp_all, log.(allα), color = "#285C93", alpha = 0.5)
#     F_n = Symbol("f", fn)
#     @eval $F_n = f
#     # save("../results/1com_$(Dnames[fn]).png", f) 
#     next!(progress)
# end

# f1
progress = Progress(N; desc="Progress running:")
temp = collect(Temp_rich .+273.15)
f1 = Figure(size = (1200, 1200));
@time for i in 1:N 
        αii = [all_ℵii[t][i] for t in 1:num_temps]
        Nα, init_in, AIC_in, temp_all, allα = try_params(αii, num_temps, 2000)
        fit_ii = curve_fit(temp_SS, temp, αii, init_in)
        ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
        scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
        lines!(ax1, Temp_rich, abs.(temp_SS(temp, fit_ii.param)), color = ("#E17542", 1), linewidth = 1)
        next!(progress)
    end 

f1
R"library(beepr); beep(sound = 4, expr = NULL)"
# save("/Users/Danica/Documents/temp_interactions/results/αii_fitted.png", f1) 

