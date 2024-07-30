
N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]

tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

T = 273.15 + 15
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
p_lv = Eff_LV_params(p=p, sol=sol);


Ci = fill(0.1, N)
prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
C = sol_LV.u[length(sol_LV.t)][1:N]

LV_Jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
bm = sol.u[length(sol.t)][1:p_lv.N]
sur = (1:p_lv.N)[bm .> 1.0e-7]

Jac_off = [sum(LV_Jac[j,i] for j in 1:N if j != i) for i in 1:N ]
diag(LV_Jac) - Jac_off

sum(diag(LV_Jac))
sum(eigen(LV_Jac).values)
sum(diag(LV_Jac) - Jac_off .> 0)



bm = sol.u[length(sol.t)][1:p_lv.N]
sur = (1:p_lv.N)[bm .> 1.0e-7]
N_s = length(sur)

u_sur = p.u[sur,:]
R_t = sol.u[length(sol.t)][N+1:N+M]
C_t = sol.u[length(sol.t)][1:N][sur]
u_tR = mapslices(x -> x .* R_t, u_sur, dims=2) # getting the actual uptake
u_t = mapslices(x -> x .* C_t, u_tR, dims=1) # getting the actual uptake

R_over = 1 .-[bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
R_over_matrix = reshape(1 .-[bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s], N_s, N_s)

l_t = p.l[sur,:,:]
ul = zeros(Float64, N_s, M)
for s in 1: N_s
    uli = zeros(Float64, M, M)
    for α in 1:M
        uli[α,:] = u_t[s, α] .* l_t[s, α, :]
    end 
    ul[s,:] = sum(uli, dims = 1)
end 
ul_over = 1 .- [bray_curtis_dissimilarity(ul[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
ul_over_matrix = transpose(reshape(1 .-[bray_curtis_dissimilarity(ul[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s], N_s, N_s))

est = ul_over - R_over
Plots.histogram(ul_over./R_over)
Plots.histogram(R_over)
Plots.histogram(ul_over)
Plots.histogram(ul_over./R_over)
# sur_ℵ = reshape(sur_ℵ, N_sur*N_sur, 1)
# Plots.scatter(est_ℵ, sur_ℵ)
Plots.histogram(est)