include("./sim_frame.jl");

N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
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

eigen(LV_jac).values

#################  

N=100
M=50
L = 0.3
### Temp params 
num_temps = 31
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)


Eff_results = zeros(Float64, num_temps, 20)
@time for i in range(0, stop = num_temps, length = num_temps+1)
    T = 273.15 + i
    all_ℵ = Float64[]; ℵii = Float64[]; ℵij = Union{Float64, Missing}[]; all_r = Float64[]; 
    all_u = Float64[]; all_m = Float64[]; RO = Union{Float64, Missing}[]; all_Eu = Float64[]; all_Em = Float64[]
    for j in 1:50 
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        ## run simulation
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        ## getting effective LV coefficients
        p_lv = Eff_LV_params(p=p, sol=sol);
        # number of species with r>0 at equilibium 
        N_sur = sum(p_lv.r .> 0)
        sur_r = p_lv.r[p_lv.r.>0]
        mean_r = mean(sur_r) 
        # mean uptake and respiration 
        m = mean(p.m[p_lv.r.>0])
        u = mean(sum(p.u, dims =2)[p_lv.r.>0])
        # mean E Tp for u and m 
        Eu, Em = mean(p.E[p_lv.r.>0, :],dims = 1)
        # Resource overlap
        u_sur = p.u[p_lv.r .> 0,:]
        # survivor interaction coefficients
        sur_ℵ = p_lv.ℵ[p_lv.r.>0, p_lv.r.>0] # interactions in the possibily surviving community
        mean_ℵ = mean(sur_ℵ) # average interaction in the possibly survived community
        mean_ℵii = mean(diag(sur_ℵ))
        # eigenvalue for jacobian 
        LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
        ev = sum(eigen(LV_jac).values .< 0)
        if N_sur > 1
            R_over = mean([cosine_similarity(u_sur[i,:], u_sur[j,:]) for i in 1:N_sur for j in 1:N_sur if i != j])
            mean_ℵij = mean([sur_ℵ[i, j]/diag(sur_ℵ)[i] for i in 1:N_sur for j in 1:N_sur if i != j])
            push!(all_ℵ, mean_ℵ); push!(ℵii, mean_ℵii); push!(ℵij,mean_ℵij); push!(all_r, mean_r);
            push!(all_u, u); push!(all_m, m); push!(RO, R_over); push!(all_Eu, Eu); push!(all_Em, Em); 
            push!(less0_ev, ev)
        else 
            push!(all_ℵ, mean_ℵ); push!(ℵii, mean_ℵii); push!(ℵij,missing); push!(all_r, mean_r);
            push!(all_u, u); push!(all_m, m); push!(RO, missing); push!(all_Eu, Eu); push!(all_Em, Em); 
            push!(less0_ev, ev)
        end
    end 
    Eff_results[Int(i+1),:] = [mean(all_ℵ), std(all_ℵ)/sqrt(length(all_ℵ)), 
        mean(ℵii), std(ℵii)/sqrt(length(ℵii)), mean(skipmissing(ℵij)), std(skipmissing(ℵij))/sqrt(length(ℵij)), 
        mean(all_r), std(all_r)/sqrt(length(all_r)), mean(all_u), std(all_u)/sqrt(length(all_u)), 
        mean(all_m), std(all_m)/sqrt(length(all_m)), mean(skipmissing(RO)), std(skipmissing(RO))/sqrt(length(RO)),
        mean(all_Eu), std(all_Eu)/sqrt(length(all_Eu)), mean(all_Em), std(all_Em)/sqrt(length(all_Em)),
        mean(less0_ev), std(less0_ev)/sqrt(length(less0_ev))]
    print(i, " °C Complete, ", "α ",mean(all_ℵ),"\n") 
end

col_names_EF = ["α", "α_err", "αii", "αii_err", "αij", "αij_err", 
                   "r", "r_err", "u", "u_err","m", "m_err", "RO", "RO_err",
                   "Eu", "Eu_err", "Em", "Em_err", "eigen", "eigen_err"];
Eff_results = DataFrame(Eff_results, col_names_EF);


Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αii", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,1], ylabel = "αij", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, Eff_results.αii, color = ("#FA8328",0.8), linewidth = 5, label = "MiCRM simulation")
band!(ax1, Temp_rich, Eff_results.αii .- Eff_results.αii_err, Eff_results.αii .+ Eff_results.αii_err, color = ("#FA8328", 0.2))
lines!(ax2, Temp_rich, Eff_results.αij, color = ("#015845", 0.8), linewidth = 5, label = "CUE Variance")
band!(ax2, Temp_rich,  Eff_results.αij .- Eff_results.αij_err, Eff_results.αij .+ Eff_results.αij_err, color = ("#015845", 0.2))
linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#FA8328",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#015845", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
# Label(f[1,1, TopLeft()], "(a)")
f

############# Distributions 
using Plots

num_temps = 7
UD_LD = []; pUD = []; pLD = []; rUDLD = []

@time for i in 1:num_temps
    T = 273.15 + 5*(i-1)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

    p_lv = Eff_LV_params(p=p, sol=sol);

    ### for getting the relative α
    # N_sur = sum(p_lv.r .> 0)
    # sur_ℵ = all_ℵ[p_lv.r .>0, p_lv.r .>0] # interactions in the possibily surviving community
    UD = [p_lv.ℵ[i, j]/ diag(p_lv.ℵ)[i] for i in 1:N for j in 1:N if i > j]
    LD = [p_lv.ℵ[i, j]/ diag(p_lv.ℵ)[i] for i in 1:N for j in 1:N if i < j]
    r_UDLD = UD./LD
    p_rUDLD = histogram(r_UDLD, title = "T = $(5*(i-1)) °C", bins = 50, color = "#6B8EDE", xlabel = "UD/LD", ylabel = "frequency", legend = false , size = (800, 800))
    p_UDLD = Plots.scatter(LD, UD, title = "T = $(5*(i-1)) °C", xlabel = "LD", 
        ylabel = "UD", color = "#6B8EDE", markerstrokewidth = 1, markerstrokecolor = "#283747", legend = false, 
        xlims = (0, 20), ylims = (0, 20),
        # xlims = (-0.05, 0), ylims = (-0.05, 0),
        size = (600, 900))
        Plots.plot!(0:20, 0:20, color ="#4F363E")
    p_UDp = histogram(UD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#015845", xlabel = "UD", ylabel = "frequency", legend = false , size = (800, 800))
    p_LDp = histogram(LD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#FA8328", xlabel = "LD", ylabel = "frequency", legend = false, size = (800, 800))
    push!(rUDLD, p_rUDLD)
    push!(UD_LD, p_UDLD)
    push!(pUD, p_UDp)
    push!(pLD, p_LDp)
end 

Plots.plot(rUDLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/relative_rUDLD.png")

Plots.plot(UD_LD..., layout = (3, 3), size = (1200, 900))
savefig("../results/relative_UDLD.png")

Plots.plot(pUD..., layout = (3, 3), size = (1200, 900))
savefig("../results/relative_pUD.png")

Plots.plot(pLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/relative_pLD.png")
