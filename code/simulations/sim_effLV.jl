include("./sim_frame.jl")

N=100
M=50
L = fill(0.3, N)

### Temp params 
ρ_t= [0.0000 0.0000]; # realistic covariance
Tr=273.15+10; Ed=3.5
T = 273.15+ 0
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
Ci = fill(0.1, N)

condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

######################### Running the MCM 
Random.seed!(123)
# Random.seed!(4)
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, niche = niche)

prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
bm_all = [sol.u[n][1:N] for n in 1:length(sol.u)]
bm = reshape(vcat(bm_all...), N, length(bm_all))'

# rbm = sol.u[length(sol.t)][N+1:N+M]
######################### Running the Effective LV
p_lv = Eff_LV_params(p=p, sol=sol);
ℵii = diag(p_lv.ℵ); mean(ℵii)
ℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i != j]; mean(ℵij)

## running LV
prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
bm_LV = reshape(vcat(sol_LV.u...), N, length(sol_LV.u))'

f = Figure(fontsize = 35, size = (2000, 900));
ax1 = Axis(f[1,1], xlabel = "Time", ylabel = "Consumer Abundance", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,2], xlabel = "Time", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
for i in 1:N 
    ### to plot 0
    lines!(ax1, sol.t[1:70], bm[1:70,i], color = ("#0758AE",0.8), linewidth = 3)
    lines!(ax2, sol_LV.t[1:50], bm_LV[1:50,i], color = ("#4F363E", 0.8), linewidth = 3)
    # ### to plot 30
    # lines!(ax1, sol.t[1:60], bm[1:60,i], color = ("#0758AE",0.8), linewidth = 3)
    # lines!(ax2, sol_LV.t[1:55], bm_LV[1:55,i], color = ("#4F363E", 0.8), linewidth = 3)
end
l1 = [LineElement(color = ("#0758AE",0.8), linestyle = nothing, linewidth = 3)]
l2 = [LineElement(color = ("#4F363E", 0.8), linestyle = nothing, linewidth = 3)]
Legend(f[1,1], [l1], tellheight = false, tellwidth = false, ["MCM"], halign = :right, valign = :top)
Legend(f[1,2], [l2], tellheight = false, tellwidth = false, ["ELV"], halign = :right, valign = :top)
Label(f[1,1, TopLeft()], "(a)")
f
save("../results/MCM_ELV_dynamics0.pdf", f) 

