include("./sim_frame.jl");
using ProgressMeter, RCall
using ColorSchemes

N=100
M=50
L = fill(0.3, N)
### Temp params 
num_temps = 38
# ρ_t= [0.0000 0.0000] #[-0.3500 -0.3500]; # realistic covariance
ρ_t= [-0.9999 -0.9999] #[-0.3500 -0.3500]; # realistic covariance
niche = fill(1.0, M, N)

Tr=273.15+10; Ed=3.5 
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)


all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); 
all_up_ℵij = Vector{Vector{Float64}}(); all_low_ℵij = Vector{Vector{Float64}}(); 
all_ℵij_sum = Vector{Vector{Float64}}(); all_D_ℵij = Vector{Vector{Float64}}();
all_ℵii_sur = Vector{Vector{Float64}}(); all_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); 
all_up_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); all_low_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); 
all_ℵij_sum_sur = Vector{Vector{Union{Float64, Missing}}}(); all_D_ℵij_sur =  Vector{Vector{Union{Float64, Missing}}}();
all_R_over = Vector{Vector{Float64}}(); all_ul_over = Vector{Vector{Float64}}(); all_Rul_over = Vector{Vector{Float64}}();
all_R_over_sur = Vector{Vector{Union{Float64, Missing}}}(); all_ul_over_sur = Vector{Vector{Union{Float64, Missing}}}()
progress = Progress(num_temps; desc="Progress running:")
for i in range(0, stop = num_temps-1, length = num_temps)
    T = 273.15 + i
    Random.seed!(6)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, niche = niche)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    sur = (1:N)[bm .> 1.0e-7]
    
    u_all = mapslices(x -> x .* x0[1:N], p.u, dims=1) # getting the actual uptake
    # R_over = 1 .- [bray_curtis_dissimilarity(u_all[i,:], u_all[j,:]) for i in 1:N for j in 1:N if j != i]
    R_over = [Jaccard_Index(u_all[i,:], u_all[j,:]) for i in 1:N for j in 1:N if j != i]
    # R_over = 1 ./(1 .+ R_eu) #; mean(R_over)
    ul = zeros(Float64, N, M)
    for s in 1: N
        uli = zeros(Float64, M, M)
        for α in 1:M
        uli[α,:] = u_all[s, α] .* p.l[s, α, :]
        end 
        ul[s,:] = sum(uli, dims = 1)
    end 
    ul_over = [Jaccard_Index(u_all[i,:], ul[j,:]) for i in 1:N for j in 1:N if j != i]
    # ul_over = 1 .- [bray_curtis_dissimilarity(u_all[i,:], ul[j,:]) for i in 1:N for j in 1:N if j != i]
    Rul_over = ul_over - R_over
    mean(Rul_over)
    sur = (1:N)[bm .> 1.0e-7]
    u_sur = p.u[sur,:]
    N_s = length(sur)
    R_t = sol.u[length(sol.t)][N+1:N+M]
    C_t = bm[bm .> 1.0e-7]

    # println(bm .> 1.0e-7)
    p_lv = Eff_LV_params(p=p, sol=sol);
    # number of species with r>0 at equilibium 
    ℵ = p_lv.ℵ
    ℵii = diag(ℵ)
    ℵij_sum = vec(sum(ℵ, dims = 2).- diag(ℵ))
    ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i != j]
    up_ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i > j]
    low_ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i < j]
    D_ℵij = up_ℵij./low_ℵij

    sur_ℵ = ℵ[sur, sur]
    ℵii_sur = diag(sur_ℵ)
    if N_s > 1
        ℵij_sum_sur = vec(sum(sur_ℵ, dims = 2).- diag(sur_ℵ))
        ℵij_sur = [sur_ℵ[i,j] for i in 1:N_s for j in 1:N_s if i != j]
        up_ℵij_sur = [sur_ℵ[i,j] for i in 1:N_s for j in 1:N_s if i > j]
        low_ℵij_sur = [sur_ℵ[i,j] for i in 1:N_s for j in 1:N_s if i < j]
        D_ℵij_sur = up_ℵij_sur./low_ℵij_sur

        u_tR = mapslices(x -> x .* R_t, u_sur, dims=2) # getting the actual uptake
        u_t = mapslices(x -> x .* C_t, u_tR, dims=1) # getting the actual uptake
        R_over_sur = 1 ./ (1 .+ [euclidean_distance(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i])
        l_t = p.l[sur,:,:]
        ul_sur = zeros(Float64, N_s, M)
        for s in 1: N_s
                uli = zeros(Float64, M, M)
                for α in 1:M
                uli[α,:] = u_t[s, α] .* l_t[s, α, :]
                end 
                ul_sur[s,:] = sum(uli, dims = 1)
        end 
        ul_over_sur = 1 ./ (1 .+ [euclidean_distance(u_t[i,:], ul_sur[j,:]) for i in 1:N_s for j in 1:N_s if j != i])
        # Rul_over_sur = ul_over_sur - R_over_sur

        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_up_ℵij, up_ℵij); push!(all_low_ℵij, low_ℵij); 
        push!(all_ℵij_sum, ℵij_sum); push!(all_D_ℵij, D_ℵij);
        push!(all_ℵii_sur, ℵii_sur); push!(all_ℵij_sur, ℵij_sur); push!(all_up_ℵij_sur, up_ℵij_sur); push!(all_low_ℵij_sur, low_ℵij_sur); 
        push!(all_ℵij_sum_sur, ℵij_sum_sur); push!(all_D_ℵij_sur, D_ℵij_sur);
        push!(all_R_over, R_over); push!(all_ul_over, ul_over); push!(all_Rul_over, Rul_over);
        push!(all_R_over_sur, R_over_sur); push!(all_ul_over_sur, ul_over_sur ); 
    else 
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_up_ℵij, up_ℵij); push!(all_low_ℵij, low_ℵij); 
        push!(all_ℵij_sum, ℵij_sum); push!(all_D_ℵij, D_ℵij);
        push!(all_ℵii_sur, ℵii_sur); push!(all_ℵij_sur, [missing]); push!(all_up_ℵij_sur, [missing]); push!(all_low_ℵij_sur, [missing]); 
        push!(all_ℵij_sum_sur, [missing]); push!(all_D_ℵij_sur, [missing]);
        push!(all_R_over, R_over); push!(all_ul_over, ul_over); push!(all_Rul_over, Rul_over);
        push!(all_R_over_sur, [missing]); push!(all_ul_over_sur, [missing]); 
    end
    next!(progress)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"


# @load "../data/1com0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur
# D = (all_ℵii = all_ℵii, all_ℵij = all_ℵij, all_up_ℵij = all_up_ℵij,all_low_ℵij = all_low_ℵij, all_ℵij_sum= all_ℵij_sum, all_D_ℵij = all_D_ℵij,
#     all_ℵii_sur = all_ℵii_sur,  all_ℵij_sur = all_ℵij_sur, all_up_ℵij_sur = all_up_ℵij_sur, all_low_ℵij_sur = all_low_ℵij_sur, all_ℵij_sum_sur = all_ℵij_sum_sur, all_D_ℵij_sur = all_D_ℵij_sur);
# Dnames = ("αii", "αij", "up_αij", "low_αij", "sum_αij", "up_low", "αii_sur", "αij_sur", "up_αij_sur", "low_αij_sur", "sum_αij_sur", "up_low_sur");

@save "../data/1com0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur all_R_over all_ul_over all_R_over_sur all_ul_over_sur

# @load "../data/1com_0.jld2" all_ℵii all_ℵij all_up_ℵij all_low_ℵij all_ℵij_sum all_D_ℵij all_ℵii_sur all_ℵij_sur all_up_ℵij_sur all_low_ℵij_sur all_ℵij_sum_sur all_D_ℵij_sur all_R_over all_ul_over all_R_over_sur all_ul_over_sur

function ellipse_points(center, axis_lengths, angle_ellipse; n_points=100)
    t = range(0, 2π, length=n_points)
    ellipse = [axis_lengths[1] .* cos.(t), axis_lengths[2] .* sin.(t)]  # scale to eigenvalues
    rotation_matrix = [cos(angle_ellipse) -sin(angle_ellipse); sin(angle_ellipse) cos(angle_ellipse)]
    rotated_ellipse = rotation_matrix * ellipse
    [rotated_ellipse[i] .+ center[i] for i in 1:2] # Center vector broadcasted across each column of the rotated ellipse
end

num_temps = 31
cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])
f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "log(|αᵢⱼ|)", ylabel = "log(|αⱼᵢ|)", xlabelsize = 50, ylabelsize = 50)
for i in 1: 31
    x = log.(abs.(all_up_ℵij[i])); y = log.(abs.(all_low_ℵij[i]))
    mean_x, mean_y = mean(x), mean(y)
    eigen_vals, eigen_vecs = eigen(cov(hcat(x, y)))
    axis_lengths = sqrt.(eigen_vals) * 2
    angle_ellipse = atan(eigen_vecs[2, 1], eigen_vecs[1, 1])
    ellipse_pts = ellipse_points([mean_x, mean_y], axis_lengths, angle_ellipse)
    scatter!(ax, log.(abs.(all_up_ℵij[i])), log.(abs.(all_low_ℵij[i])), color = cs[i], markersize = 5, alpha = 0.2)
    lines!(ax, ellipse_pts[1], ellipse_pts[2], color=cs[i], linewidth=2)
end 
min = minimum(log.(abs.(vcat(all_up_ℵij..., all_low_ℵij...))))
max = maximum(log.(abs.(vcat(all_up_ℵij..., all_low_ℵij...))))
lines!(ax, [min, max], [min, max], linestyle = :dash, color = ("#4F363E", 0.9), linewidth = 2)
Colorbar(f[1,2], colorrange = [0, num_temps], colormap = cs, label = "Temperature")
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/αijαji-1.pdf", f) 

# i = 11
# f = Figure(fontsize = 35, size = (1200, 900));
# ax = Axis(f[1,1], xlabel = "log(|αij|)", ylabel = "log(|αji|)", xlabelsize = 50, ylabelsize = 50)
# scatter!(ax, log.(abs.(all_up_ℵij[i])), log.(abs.(all_low_ℵij[i])), color = cs[i], markersize = 15, alpha = 0.1)
# Colorbar(f[1,2], colorrange = [0, num_temps], colormap = cs, label = "Temperature")
# f

Temp_rich = range(0, num_temps-1, length = num_temps)
k = 0.0000862
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)
temp_αii = vcat([repeat([Temp_rich[t]], length(all_ℵii[t])) for t in 1:num_temps]...)

allαii = log.(abs.(vcat(all_ℵii...)))

mean_R_over = [mean(all_R_over[t]) for t in 1:num_temps]
mean_R_err = [std(all_R_over[t])/sqrt(length(all_R_over[t])) for t in 1:num_temps]
mean_ul = [mean(all_ul_over[t]) for t in 1:num_temps]
mean_ul_err = [std(all_ul_over[t])/sqrt(length(all_ul_over[t])) for t in 1:num_temps]
mean_est = [mean(all_Rul_over[t]) for t in 1:num_temps]
mean_est_err = [std(all_Rul_over[t])/sqrt(length(all_Rul_over[t])) for t in 1:num_temps]

Temp_rich = range(0, 30, length = 31)
# temp_R_over = vcat([repeat([Temp_rich[t]], length(all_R_over[t])) for t in 1:31]...)
# temp_ul_over = vcat([repeat([Temp_rich[t]], length(all_ul_over[t])) for t in 1:31]...)
f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Competition & Facilitation", ygridvisible = false, xgridvisible = false, xlabelsize = 50, ylabelsize = 50)
ax2 = Axis(f[1,1], ylabel = "Pairwise Interaction", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, mean_R_over, linewidth = 5, color = ("#285C93", 0.7), label = "ƒo")
band!(ax1, Temp_rich, mean_R_over .- mean_R_err , mean_R_over.+ mean_R_err, color = ("#376298", 0.3))
lines!(ax1, Temp_rich, mean_ul, linewidth = 5, color = ("#E17542", 0.7), label = "ƒc")
band!(ax1, Temp_rich, mean_ul .- mean_ul_err , mean_ul .+ mean_ul_err , color = ("#E17542", 0.3))
# lines!(ax2, Temp_rich, mean_est, color = ("#4F363E",0.8), linewidth = 5, label = "")
# band!(ax2, Temp_rich, mean_estα .- mean_est_err , mean_est .+ mean_est_err, color = ("#4F363E", 0.3))
# boxplot!(ax, temp_R_over, vcat(all_R_over[1:31]...), color = ("#285C93", 0.5), label = "ƒo")
# boxplot!(ax, temp_ul_over, vcat(all_ul_over[1:31]...), color = ("#E17542", 0.5), label = "ƒc")
# axislegend(position = :rt)
l1 = [LineElement(color = ("#285C93", 0.7), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#E17542", 0.7), linestyle = nothing, linewidth = 5)]
l3 = [LineElement(color = ("#4F363E", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2, l3], tellheight = false, tellwidth = false, [ "ƒo", "ƒc", "ƒc-ƒo"], halign = :right, valign = :top)
Label(f[1,1, TopLeft()], "(c)")
f
save("../results/comp_facT0.pdf", f) 

cscheme = ColorScheme(range(colorant"#376298",colorant"#ECDFCB", length = 16))
cscheme1 = ColorScheme(range(colorant"#ECDFCB",colorant"#9A2B1A", length = 16))
cs = vcat(cscheme[1:16], cscheme1[2:16])
f = Figure(fontsize = 35, size = (1200, 1200));
ax = Axis(f[1,1], xlabel = "Competition", ylabel = "Facilitation", xlabelsize = 50, ylabelsize = 50)
for i in 1: 31
    scatter!(ax, all_R_over[i], all_ul_over[i], color = cs[i], markersize = 5, alpha = 0.2)
end 
# scatter!(ax, real.(vcat(all_circ...)), vcat(all_leadH...), color = ("#285C93"), label = "", markersize = 10, alpha = 0.3) 
# lines!(ax, [ minimum(real.(vcat(all_circ...))), maximum(real.(vcat(all_circ...)))], [0,0], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
lines!(ax, [minimum(vcat(all_R_over...)[vcat(all_R_over...).<1]),maximum(vcat(all_R_over...)[vcat(all_R_over...).<1])], [minimum(vcat(all_R_over...)[vcat(all_R_over...).<1]), maximum(vcat(all_R_over...)[vcat(all_R_over...).<1])], linestyle = :dash, color = ("#4F363E", 1), linewidth = 3)
Colorbar(f[1,2], colorrange = [0, 30], colormap = cs, label = "Temperature")
Label(f[1,1, TopLeft()], "(b)")
f
save("../results/comp_fac-1.pdf", f) 

temp_ℵij = vcat([repeat([Temp_rich[t]], length(all_ℵij[t])) for t in 1:31]...)

all_ℵij
f = Figure(fontsize = 35,size = (1200, 1200));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Competition & Facilitation", ygridvisible = false, xgridvisible = false, xlabelsize = 50, ylabelsize = 50)
for i in 1: 31
    # scatter!(ax, all_R_over[i], all_ul_over[i], color = cs[i], markersize = 5, alpha = 0.1)
    hist!(ax, all_ℵij[i], color = cs[i], label = "Competition", bins = 100)

end 
Label(f[1,1, TopLeft()], "(a)")
f





########### Fitting community α #################
include("./fitting.jl");

# f1
progress = Progress(N; desc="Progress running:")
temp = collect(Temp_rich .+273.15)
f1 = Figure(size = (1200, 1200));
fitted = zeros(Float64, N, 6)
@time for i in 1:N 
        next!(progress)
        αii = [all_ℵii[t][i] for t in 1:num_temps]
        Nα, init_in, AIC_in, temp_all, allα = try_params(αii, num_temps, 2000)
        fit_ii = curve_fit(temp_SS, temp_all, allα, init_in)
        r_square = calculate_r2(fit_ii, temp_all, allα)
        params = fit_ii.param
        ## calculate_r2
        pred = abs.(temp_SS(temp, params))
        ss_res = sum((allα .- pred).^2)
        ss_tot = sum((allα .- mean(allα)).^2)
        r_square = 1 - ss_res / ss_tot
        ## calculate_AIC
        aic_value = N * log(ss_res / N) + 2 * 4
        ## store data
        fitted[Int(i),:] = vcat(params, aic_value, r_square)
        ## plotting data
        ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
        scatter!(ax1, Temp_rich, abs.(αii), color = "#285C93", alpha = 0.5)
        lines!(ax1, Temp_rich, pred, color = ("#E17542", 1), linewidth = 1)
    end 
R"library(beepr); beep(sound = 4, expr = NULL)"

df_names = ["B0","E","Th","Ed","AIC","r2"]
fitted = DataFrame(fitted, df_names);
# CSV.write("../results/αii_fitted_0.csv", fitted, writeheader=false)
CSV.write("../results/αii_fitted-1.csv", fitted, writeheader=false)

f1
# save("/Users/Danica/Documents/temp_interactions/results/αii_fitted_0.pdf", f1) 
# save("/Users/Danica/Documents/temp_interactions/results/αij_sum_fitted_0.pdf", f1) 
# fitted = CSV.read("../results/αii_fitted-1.csv", DataFrame, header=false)

fitted = fitted[fitted.E .> 0, :]
fitted = fitted[fitted.B0 .> 0, :]

mean_Bii = mean(log.(fitted.B0))
var_Bii = var(log.(fitted.B0))/abs(mean(log.(fitted.B0)))
mean_Eii = mean(fitted.E)
var_Eii = var(fitted.E)/abs(mean(fitted.E))
cor_ii = cor(log.(fitted.B0), fitted.E)

Nα, B_m, E_m, T_m, Ed,temp_all, allα = get_init_param(all_ℵii, num_temps)
E_p = mean(fitted.E[fitted.E .> 0])
T_p = mean(fitted.Th[fitted.Th .> 0])
init_in = [exp(B_m), E_p, T_p, Ed]
fit_ii = curve_fit(temp_SS, temp_all, allα, init_in)
## calculate_r2
r_square = calculate_r2(fit_ii, temp_all, allα)
pred = abs.(temp_SS(temp, fit_ii.param))
E = fit_ii.param[2]

f = Figure(fontsize = 35,size = (1200, 1200));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "log(|αij|)", ygridvisible = false, xgridvisible = false, xlabelsize = 50, ylabelsize = 50)
scatter!(ax, temp_all .- 273.15, log.(abs.(allα)), color = "#285C93", alpha = 0.5)
lines!(ax, Temp_rich, log.(pred), color = ("#E17542", 1), linewidth = 1)
text!(10, -4.5, text = "E = $(round(E, digits = 3))", align = (:center, :center), fontsize = 35)
f
# save("/Users/Danica/Documents/temp_interactions/results/αii_fitted_all-1.pdf", f) 
