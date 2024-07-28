include("./sim_frame.jl");
using Plots

N=7
M=5
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

############# Distributions 

num_temps = 7
UD_LD = []; pUD = []; pLD = []; rUDLD = []; 
rela_UD_LD = []; rela_pUD = []; rela_pLD = []; rela_rUDLD = [];
pr =[]; pr_sur = []; pαii = []; pαii_sur =[]; pαij = []; pαij_sur = []

@time for i in 1:num_temps
    T = 273.15 + 5*(i-1)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    p_lv = Eff_LV_params(p=p, sol=sol);

    #### Calculating Jacobian
    LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
    dia_Jac = diag(LV_jac)
    # leading = maximum(real(eigen(LV_jac).values))

    ### for getting the relative α
    N_sur = sum(p_lv.r .> 0)
    αii = histogram(log.(abs.(diag(p_lv.ℵ))), title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#5B8867", xlabel = "αii", ylabel = "frequency", legend = false , size = (800, 800))
    αii_sur = histogram(log.(abs.(diag(p_lv.ℵ)[p_lv.r.>0])), title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#EF8F8C", xlabel = "αii_sur", ylabel = "frequency", legend = false , size = (800, 800))
    αij = histogram(log.(abs.([p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i != j])), title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#5B8867", xlabel = "αij", ylabel = "frequency", legend = false , size = (800, 800))
    αij_sur = histogram(log.(abs.([p_lv.ℵ[p_lv.r.>0, p_lv.r.>0][i, j] for i in 1:N_sur for j in 1:N_sur if i != j])), title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#EF8F8C", xlabel = "αij_sur", ylabel = "frequency", legend = false , size = (800, 800))
    r = histogram(log.(abs.(p_lv.r)), title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#5B8867", xlabel = "r", ylabel = "frequency", legend = false , size = (800, 800))
    r_sur = histogram(log.(p_lv.r[p_lv.r.>0]), title = "T = $(5*(i-1)) °C", bins = 50, 
                color = "#EF8F8C", xlabel = "r_sur", ylabel = "frequency", legend = false , size = (800, 800))
    ### absolute alphas ###
    UD = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i > j]
    LD = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i < j]
    r_UDLD = UD./LD
    p_rUDLD = histogram(r_UDLD, title = "T = $(5*(i-1)) °C", bins = 50, color = "#6B8EDE", xlabel = "UD/LD", ylabel = "frequency", legend = false , size = (800, 800))
    p_UDLD = Plots.scatter(LD, UD, title = "T = $(5*(i-1)) °C", xlabel = "LD", 
        ylabel = "UD", color = "#6B8EDE", markerstrokewidth = 1, markerstrokecolor = "#283747", legend = false, 
        # xlims = (0, 20), ylims = (0, 20),
        xlims = (-0.05, 0), ylims = (-0.05, 0),
        size = (600, 900))
        Plots.plot!(-0.05:0, -0.05:0, color ="#4F363E")
    p_UDp = histogram(UD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#015845", xlabel = "UD", ylabel = "frequency", legend = false , size = (800, 800))
    p_LDp = histogram(LD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#FA8328", xlabel = "LD", ylabel = "frequency", legend = false, size = (800, 800))

    ###relative alphas ###
    rela_UD = [p_lv.ℵ[i, j]/ diag(p_lv.ℵ)[i] for i in 1:N for j in 1:N if i > j]
    rela_LD = [p_lv.ℵ[i, j]/ diag(p_lv.ℵ)[i] for i in 1:N for j in 1:N if i < j]
    rela_r_UDLD = rela_UD./rela_LD
    rela_p_rUDLD = histogram(rela_r_UDLD, title = "T = $(5*(i-1)) °C", bins = 50, color = "#6B8EDE", xlabel = "UD/LD", ylabel = "frequency", legend = false , size = (800, 800))
    rela_p_UDLD = Plots.scatter(rela_LD, rela_UD, title = "T = $(5*(i-1)) °C", xlabel = "LD", 
        ylabel = "UD", color = "#6B8EDE", markerstrokewidth = 1, markerstrokecolor = "#283747", legend = false, 
        xlims = (0, 20), ylims = (0, 20),
        # xlims = (-0.05, 0), ylims = (-0.05, 0),
        size = (600, 900))
        Plots.plot!(0:20, 0:20, color ="#4F363E")
    rela_p_UDp = histogram(rela_UD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#015845", xlabel = "UD", ylabel = "frequency", legend = false , size = (800, 800))
    rela_p_LDp = histogram(rela_LD , title = "T = $(5*(i-1)) °C", bins = 50, color = "#FA8328", xlabel = "LD", ylabel = "frequency", legend = false, size = (800, 800))
    
    ### Collecting plots ###
    push!(rUDLD, p_rUDLD); push!(UD_LD, p_UDLD); push!(pUD, p_UDp); push!(pLD, p_LDp)
    push!(rela_rUDLD, rela_p_rUDLD); push!(rela_UD_LD, rela_p_UDLD); push!(rela_pUD, rela_p_UDp); push!(rela_pLD, rela_p_LDp)
    push!(pr, r); push!(pr_sur, r_sur); push!(pαii, αii); push!(pαii_sur, αii_sur); push!(pαij, αij); push!(pαij_sur, αij_sur)
end 

Plots.plot(pr..., layout = (3, 3), size = (1200, 900))
savefig("../results/pr.png")
Plots.plot(pr_sur..., layout = (3, 3), size = (1200, 900))
savefig("../results/pr_sur.png")
Plots.plot(pαii..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαii.png")
Plots.plot(pαii_sur..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαii_sur.png")
Plots.plot(pαij..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαij.png")
Plots.plot(pαij_sur..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαij_sur.png")

Plots.plot(rUDLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/rUDLD.png")
Plots.plot(UD_LD..., layout = (3, 3), size = (1200, 900))
savefig("../results/UDLD.png")
Plots.plot(pUD..., layout = (3, 3), size = (1200, 900))
savefig("../results/pUD.png")
Plots.plot(pLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/pLD.png")

Plots.plot(rela_rUDLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/rela_rUDLD.png")
Plots.plot(rela_UD_LD..., layout = (3, 3), size = (1200, 900))
savefig("../results/rela_UDLD.png")
Plots.plot(rela_pUD..., layout = (3, 3), size = (1200, 900))
savefig("../results/rela_pUD.png")
Plots.plot(rela_pLD..., layout = (3, 3), size = (1200, 900))
savefig("../results/rela_pLD.png")
