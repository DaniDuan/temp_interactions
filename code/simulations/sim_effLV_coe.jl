include("./sim_frame.jl");

N=100
M=50
L = 0.3
### Temp params 
num_temps = 31
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
iter = 50

Eff_results = zeros(Float64, num_temps, 27)
@time for i in range(0, stop = 30, length = 31)
    T = 273.15 + i
    all_ℵ = Float64[]; ℵii = Float64[]; ℵij = Union{Float64, Missing}[]; all_r = Float64[]; 
    all_leading = Float64[]; stability = Float64[]; all_diag = Float64[];radi = Float64[]; diag_dominance = Float64[];
    all_u = Float64[]; all_m = Float64[]; RO = Union{Float64, Missing}[]; all_Eu = Float64[]; all_Em = Float64[]
    for j in 1:iter 
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
        jac_eigen = eigen(LV_jac).values
        leading = maximum(real.(jac_eigen))
        jac_diag = diag(LV_jac)
        jac_off = [sum(LV_jac[j, i] for j in 1:N if j != i) for i in 1:N ]
        diag_dom = sum(abs.(jac_diag) - abs.(jac_off) .> 0)/N
        if N_sur > 1
            R_over = mean([cosine_similarity(u_sur[i,:], u_sur[j,:]) for i in 1:N_sur for j in 1:N_sur if i != j])
            mean_ℵij = mean([sur_ℵ[i, j]/diag(sur_ℵ)[i] for i in 1:N_sur for j in 1:N_sur if i != j])
            push!(all_ℵ, mean_ℵ); push!(ℵii, mean_ℵii); push!(ℵij,mean_ℵij); push!(all_r, mean_r);
            push!(all_u, u); push!(all_m, m); push!(RO, R_over); push!(all_Eu, Eu); push!(all_Em, Em); 
            push!(all_leading, leading); append!(all_diag, jac_diag); append!(radi, jac_off); push!(diag_dominance, diag_dom)
        else 
            push!(all_ℵ, mean_ℵ); push!(ℵii, mean_ℵii); push!(ℵij,missing); push!(all_r, mean_r);
            push!(all_u, u); push!(all_m, m); push!(RO, missing); push!(all_Eu, Eu); push!(all_Em, Em); 
            push!(all_leading, leading); append!(all_diag, jac_diag); append!(radi, jac_off); push!(diag_dominance, diag_dom)
        end
    end 
    Eff_results[Int(i+1),:] = [mean(all_ℵ), std(all_ℵ)/sqrt(length(all_ℵ)), 
        mean(ℵii), std(ℵii)/sqrt(length(ℵii)), mean(skipmissing(ℵij)), std(skipmissing(ℵij))/sqrt(length(ℵij)), 
        mean(all_r), std(all_r)/sqrt(length(all_r)), mean(all_u), std(all_u)/sqrt(length(all_u)), 
        mean(all_m), std(all_m)/sqrt(length(all_m)), mean(skipmissing(RO)), std(skipmissing(RO))/sqrt(length(RO)),
        mean(all_Eu), std(all_Eu)/sqrt(length(all_Eu)), mean(all_Em), std(all_Em)/sqrt(length(all_Em)),
        mean(all_leading), std(all_leading)/sqrt(length(all_leading)), sum(all_leading .< 0)/iter, 
        mean(diag_dominance), std(diag_dominance)/sqrt(length(diag_dominance)), 
        mean(all_diag), std(all_diag)/sqrt(length(all_diag)), mean(radi), std(radi)/sqrt(length(radi))]
    print(i, " °C Complete, ", "α ",mean(all_ℵ),"\n") 
end

col_names_EF = ["α", "α_err", "αii", "αii_err", "αij", "αij_err", 
                   "r", "r_err", "u", "u_err","m", "m_err", "RO", "RO_err",
                   "Eu", "Eu_err", "Em", "Em_err", "eigen", "eigen_err", "stability",
                   "diag_dom", "diag_dom_err",
                   "Jac_diag", "Jac_diag_err", "radius", "radius_err"];
Eff_results = DataFrame(Eff_results, col_names_EF);

# CSV.write("../data/Eff_results.csv", Eff_results, writeheader=false)

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

