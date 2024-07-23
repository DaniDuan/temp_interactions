include("../sim_frame.jl")

################### ϵ vs. R* ###############################
N=100
M=50
l_α = 0.3
### Temp params 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+9; Ed=3.5 
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 31

ϵ_sur = Float64[]; ϵ_ext = Float64[]; R_sur = Float64[]; R_ext = Float64[]
T = 273.15 + 10
## generate params
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=l_α, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
## Calc CUE
ϵ = (p.u * x0[N+1:N+M] .* (1-l_α) .- p.m) ./ (p.u * x0[N+1:N+M])
Rs = p.m ./(p.u * x0[N+1:N+M] .* (1-l_α))
## run simulation
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N] 
ϵ_sur= ϵ[bm.>1e-7]; ϵ_ext= ϵ[bm.<=1e-7]; R_sur= Rs[bm.>1e-7]; R_ext= Rs[bm.<=1e-7]

# plots
f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "CUE", ylabel = "R*", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, ϵ_ext, R_ext, color = ("#4F363E", 0.4), markersize = 25, label = "Extinct")
scatter!(ax, ϵ_sur, R_sur, color = ("#EF8F8C",1), marker = :star4, markersize = 25, label = "Survivor")
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(b)")
f

save("../result/CUE_R.png", f) 
