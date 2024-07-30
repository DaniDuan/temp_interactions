include("./sim_frame.jl");
using Plots

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
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

############# Distributions 

num_temps = 7
UD_LD = []; pUD = []; pLD = []; rUDLD = []; 
rela_UD_LD = []; rela_pUD = []; rela_pLD = []; rela_rUDLD = []; 
pαii = []; pαii_sur =[]; pαij = []; pαij_sur = []; prela_αij = []; prela_αij_sur = [];
pp_RO = []; pp_ul = []; pest_α = [];

@time for i in 1:num_temps
    T = 273.15 + 5*(i-1)
    R_over_v = Float64[]; ul_over_v = Float64[]; est = Float64[];
    UD = Float64[]; LD = Float64[]; r_UDLD = Float64[];
    rela_UD = Float64[]; rela_LD = Float64[]; rela_r_UDLD = Float64[];
    αii_v = Float64[]; αii_v_sur = Float64[]; αij_v = Float64[]; αij_v_sur = Float64[]; 
    rela_αij_v = Float64[]; rela_αij_v_sur = Float64[];
    for j in 1: 50
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        ## run simulation
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        p_lv = Eff_LV_params(p=p, sol=sol);
        #### resource overlap and cross-feeding 
        bm = sol.u[length(sol.t)][1:p_lv.N]
        sur = (1:p_lv.N)[bm .> 1.0e-7]
        N_s = length(sur)
        N_sur = sum(p_lv.r .> 0)

        u_sur = p.u[sur,:]
        R_t = sol.u[length(sol.t)][N+1:N+M]
        C_t = sol.u[length(sol.t)][1:N][sur]
        u_tR = mapslices(x -> x .* R_t, u_sur, dims=2) # getting the actual uptake
        u_t = mapslices(x -> x .* C_t, u_tR, dims=1) # getting the actual uptake
        
        R_over = 1 .-[bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
        # R_over = 1 ./[euclidean_distance(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
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
        # ul_over = 1 ./ [euclidean_distance(ul[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
        #     ul_over_matrix = transpose(reshape(1 .-[bray_curtis_dissimilarity(ul[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s], N_s, N_s))
        
        ################ results collection #################
        # resource overlap vs. cross-feeding 
        append!(R_over_v, R_over)
        append!(ul_over_v, ul_over)
        append!(est, ul_over - R_over)
        # #### Calculating Jacobian
        # LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
        # dia_Jac = diag(LV_jac)
        # # leading = maximum(real(eigen(LV_jac).values))

        # interaction matrix distributions
        append!(αii_v, log.(abs.(diag(p_lv.ℵ))))
        append!(αii_v_sur, log.(abs.(diag(p_lv.ℵ)[sur])))
        append!(αij_v, log.(abs.([p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i != j])))
        append!(αij_v_sur, log.(abs.([p_lv.ℵ[sur, sur][i, j] for i in 1:N_s for j in 1:N_s if i != j])))
        append!(rela_αij_v, log.(abs.([p_lv.ℵ[i, j]/p_lv.ℵ[i, i] for i in 1:N for j in 1:N if i != j])))
        append!(rela_αij_v_sur, log.(abs.([p_lv.ℵ[sur, sur][i, j]/p_lv.ℵ[sur, sur][i, i] for i in 1:N_s for j in 1:N_s if i != j])))
        # original upper and lower diagonal values
        append!(UD, [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i > j])
        append!(LD, [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i < j])
        append!(r_UDLD, [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i > j]./[p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i < j])
        # upper and lower diagonal values relative to diagonal
        append!(rela_UD, [p_lv.ℵ[i, j]/ p_lv.ℵ[i, i] for i in 1:N for j in 1:N if i > j])
        append!(rela_LD, [p_lv.ℵ[i, j]/ p_lv.ℵ[i, i] for i in 1:N for j in 1:N if i < j])
        append!(rela_r_UDLD, [p_lv.ℵ[i, j]/ p_lv.ℵ[i, i] for i in 1:N for j in 1:N if i > j]./[p_lv.ℵ[i, j]/ p_lv.ℵ[i, i] for i in 1:N for j in 1:N if i < j])
    end 
    ### resource overlap vs. crossfeeding ###
    p_RO = histogram(R_over_v , title = "T = $(5*(i-1)) °C", bins = 50, color = "#015845", xlabel = "Resource Overlap", ylabel = "frequency", legend = false , size = (800, 800))
    p_ul = histogram(ul_over_v , title = "T = $(5*(i-1)) °C", bins = 50, color = "#FA8328", xlabel = "Cross-feeding", ylabel = "frequency", legend = false, size = (800, 800))
    est_α =  histogram(est, title = "T = $(5*(i-1)) °C", bins = 50, color = "#EF8F8C", xlabel = "C-R", ylabel = "frequency", legend = false , size = (800, 800))

    ### α ### 
    αii = histogram(αii_v, title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#5B8867", xlabel = "αii", ylabel = "frequency", legend = false , size = (800, 800))
    αii_sur = histogram(αii_v_sur, title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#EF8F8C", xlabel = "αii_sur", ylabel = "frequency", legend = false , size = (800, 800))
    αij = histogram(αij_v, title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#5B8867", xlabel = "αij", ylabel = "frequency", legend = false , size = (800, 800))
    αij_sur = histogram(αij_v_sur, title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#EF8F8C", xlabel = "αij_sur", ylabel = "frequency", legend = false , size = (800, 800))
    rela_αij = histogram(rela_αij_v, title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#5B8867", xlabel = "rela_αij", ylabel = "frequency", legend = false , size = (800, 800))
    rela_αij_sur = histogram(rela_αij_v_sur, title = "T = $(5*(i-1)) °C", bins = 50, 
            color = "#EF8F8C", xlabel = "rela_αij_sur", ylabel = "frequency", legend = false , size = (800, 800))

    ### upper and lower diagonal ###
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
    push!(pαii, αii); push!(pαii_sur, αii_sur); push!(pαij, αij); push!(pαij_sur, αij_sur)
    push!(prela_αij, rela_αij); push!(prela_αij_sur, rela_αij_sur);
    push!(pp_RO, p_RO); push!(pp_ul, p_ul); push!(pest_α, est_α)

    print(5*(i-1), " °C Complete, \n") 
end 

Plots.plot(pαii..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαii.png")
Plots.plot(pαii_sur..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαii_sur.png")
Plots.plot(pαij..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαij.png")
Plots.plot(pαij_sur..., layout = (3, 3), size = (1200, 900))
savefig("../results/pαij_sur.png")
Plots.plot(prela_αij..., layout = (3, 3), size = (1200, 900))
savefig("../results/prela_αij.png")
Plots.plot(prela_αij_sur..., layout = (3, 3), size = (1200, 900))
savefig("../results/prela_αij_sur.png")

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

Plots.plot(pp_RO..., layout = (3, 3), size = (1200, 900))
savefig("../results/p_RO.png")
Plots.plot(pp_ul..., layout = (3, 3), size = (1200, 900))
savefig("../results/p_ul.png")
Plots.plot(pest_α..., layout = (3, 3), size = (1200, 900))
savefig("../results/est_α.png")
