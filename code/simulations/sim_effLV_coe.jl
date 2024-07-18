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

UD = [all_ℵ[i, j] for i in 1:N for j in 1:N if i > j]
LD = [all_ℵ[i, j] for i in 1:N for j in 1:N if i < j]
plot(LD, UD)
### Removing species with growth rate lower than 0
# mean_ℵ = mean([all_ℵ[i, j] for i in 1:N for j in 1:N if i != j])/mean(diag(all_ℵ))  # the average overall interaction strength
N_sur = sum(all_r .> 0)
sur_ℵ = all_ℵ[all_r.>0, all_r.>0] # interactions in the possibily surviving community
mean_ℵ = mean(sur_ℵ) # average interaction in the possibly survived community
mean_ℵii = mean(diag(sur_ℵ))
mean_ℵij = mean([sur_ℵ[i, j] for i in 1:N_sur for j in 1:N_sur if i != j])
# mean_ℵ = mean([sur_ℵ[i, j] for i in 1:N_sur for j in 1:N_sur if i != j])/mean(diag(sur_ℵ)) # average interaction in the possibly survived community
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
ρ_t= [-0.3500, -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)


Eff_results = zeros(Float64, num_temps, 12)
@time for i in range(0, stop = 30, length = 31)
    T = 273.15 + i
    all_ℵ = Float64[]; ℵii = Float64[]; ℵij = Float64[]; all_r = Float64[]; 
    all_u = Float64[]; all_m = Float64[]
    for j in 1:50 
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        ## run simulation
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        
        p_lv = Eff_LV_params(p=p, sol=sol);
        N_sur = sum(p_lv.r .> 0)
        sur_ℵ = p_lv.ℵ[p_lv.r.>0, p_lv.r.>0] # interactions in the possibily surviving community
        mean_ℵ = mean(sur_ℵ) # average interaction in the possibly survived community
        mean_ℵii = mean(diag(sur_ℵ))
        mean_ℵij = mean([sur_ℵ[i, j] for i in 1:N_sur for j in 1:N_sur if i != j])
        sur_r = p_lv.r[p_lv.r.>0]
        mean_r = mean(sur_r) 
        m = mean(p.m[p_lv.r.>0])
        u = mean(sum(p.u, dims =2)[p_lv.r.>0])
        push!(all_ℵ, mean_ℵ); push!(ℵii, mean_ℵii); push!(ℵij,mean_ℵij); push!(all_r, mean_r);
        push!(all_u, u); push!(all_m, m)
    end 
    Eff_results[Int(i+1),:] = [mean(all_ℵ), std(all_ℵ)/sqrt(length(all_ℵ)), 
        mean(ℵii), std(ℵii)/sqrt(length(ℵii)), mean(ℵij), std(ℵij)/sqrt(length(ℵij)), 
        mean(all_r), std(all_r)/sqrt(length(all_r)), mean(all_u), std(all_u)/sqrt(length(all_u)), 
        mean(all_m), std(all_m)/sqrt(length(all_m))]
    print(i, " °C Complete, ", "α ",mean(all_ℵ),"\n") 
end 

Temp_rich = range(0, 30, 31)

f = Figure(fontsize = 30, resolution = (1200, 800));
ax = Axis(f[1,1], xlabel = "T", ylabel = "α")
scatter!(ax, Temp_rich, all_ℵ, color = ("#6B8EDE",0.8), markersize = 15)
f

############# Distributions 
using Plots

num_temps = 7
UD_LD = []; pUD = []; pLD = []

@time for i in 1:num_temps
    T = 273.15 + 5*(i-1)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

    p_lv = Eff_LV_params(p=p, sol=sol);
    UD = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i > j]
    LD = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i < j]
    p_UDLD = Plots.scatter(LD, UD, title = "T = $(5*(i-1)) °C", xlabel = "LD", 
        ylabel = "UD", color = "#6B8EDE", markerstrokewidth = 1, markerstrokecolor = "#283747", legend = false, 
        xlims = (-0.05, 0), ylims = (-0.05, 0),
        size = (800, 600))
    p_UDp = histogram(UD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#015845", xlabel = "UD", ylabel = "frequency", legend = false , size = (800, 600))
    p_LDp = histogram(LD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#FA8328", xlabel = "LD", ylabel = "frequency", legend = false, size = (800, 600))
    push!(UD_LD, p_UDLD)
    push!(pUD, p_UDp)
    push!(pLD, p_LDp)
end 

Plots.plot(UD_LD..., layout = (3, 3), size = (1200, 900))
savefig("../results/UDLD.png")

Plots.plot(pUD..., layout = (3, 3), size = (1200, 900))
savefig("../results/pUD.png")

Plots.plot(pLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/pLD.png")
