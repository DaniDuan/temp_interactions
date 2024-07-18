include("./sim_frame.jl");

N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t=[-0.1384 -0.1384]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)


################# Template for eff_LV calc ##################
T = 273.15+10
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
## run simulation
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

@time p_lv = Eff_LV_params(p=p, sol=sol);
all_ℵ = p_lv.ℵ
all_r = p_lv.r

αi = sum(all_ℵ, dims = 2) # sum effect of all species on i

### Removing species with growth rate lower than 0
# mean_ℵ = mean([all_ℵ[i, j] for i in 1:N for j in 1:N if i != j])/mean(diag(all_ℵ))  # the average overall interaction strength
N_sur = sum(all_r .> 0)
sur_ℵ = all_ℵ[all_r.>0, all_r.>0] # interactions in the possibily surviving community
mean_ℵ = mean([sur_ℵ[i, j] for i in 1:N_sur for j in 1:N_sur if i != j])/mean(diag(sur_ℵ)) # average interaction in the possibly survived community
sur_αi = sum(sur_ℵ, dims = 2) # sum effect of all species on i
sur_r = all_r[all_r.>0]
mean_r = mean(sur_r) 
m = p.m[all_r.>0]
u = sum(p.u, dims =2)[all_r.>0]

bm = sol.u[length(sol.t)][1:N]
pred = sum(bm.>1e-7)

#################  

N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t=[-0.1384 -0.1384]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

T = 273.15+10
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
## run simulation
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

@time p_lv = Eff_LV_params(p=p, sol=sol);
all_ℵ = p_lv.ℵ
all_r = p_lv.r



dead = Float64[]
@time for i in range(0, stop = 30, length = 31)
    T = 273.15 + i
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    
    p_lv = Eff_LV_params(p=p, sol=sol);
    all_ℵ = p_lv.ℵ
    all_r = p_lv.r
    push!(dead, sum(all_r .< 0))
end 

Temp_rich = range(0, 30, 7)

plot(Temp_rich, 100 .- dead)

f = Figure(fontsize = 30, resolution = (1200, 800));
ax = Axis(f[1,1], xlabel = "Relative Richness (LV)", ylabel = "Relative Richness (MiCRM)")
scatter!(ax, all_pred./maximum(all_pred), all./maximum(all), color = "#601210", markersize = 15)
f

