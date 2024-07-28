################# Template for eff_LV calc ##################
N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(0.1, M))
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

UD = [all_ℵ[i, j] for i in 1:N for j in 1:N if i > j]
LD = [all_ℵ[i, j] for i in 1:N for j in 1:N if i < j]
### Removing species with growth rate lower than 0
N_sur = sum(all_r .> 0)

# Calculate resource
u_sur = p.u[all_r .> 0,:]
R_over = mean([cosine_similarity(u_sur[i,:], u_sur[j,:]) for i in 1:N_sur for j in 1:N_sur if i != j])

sur_ℵ = all_ℵ[all_r.>0, all_r.>0] # interactions in the possibily surviving community
mean_ℵ = mean(sur_ℵ) # average interaction in the possibly survived community
mean_ℵii = mean(diag(sur_ℵ))
mean_ℵij = mean([sur_ℵ[i, j] for i in 1:N_sur for j in 1:N_sur if i != j])

relative_ℵ = mean([sur_ℵ[i, j]/diag(sur_ℵ)[i] for i in 1:N_sur for j in 1:N_sur if i != j]) # average relative interaction 

sur_r = all_r[all_r.>0]
mean_r = mean(sur_r) 
m = p.m[all_r.>0]
u = sum(p.u, dims =2)[all_r.>0]

bm = sol.u[length(sol.t)][1:N]
pred = sum(bm.>1e-7)

#### Calculating Jacobian
LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)

maximum(eigen(LV_jac).values)

