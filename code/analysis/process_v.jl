include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob

N=100
M=50
L = fill(0.3, N)
var_B = range(0, 0.3, 50)
niche = fill(1.0, M, N)

###################################
# Generate MiCRM parameters
tspan = (0.0, 2.5e10)
x0 = vcat(fill(0.1, N), fill(1, M)) 
condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

path = glob("v_*", "../data/v_notlog/")
progress = Progress(length(path)*50; desc="Progress running:")
# idx = collect(CartesianIndices(zeros(Float64, N, N)))
# ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
all_leading_collect = Vector{Vector{ComplexF64}}(); all_rich_collect = Vector{Vector{Float64}}(); 
all_r_sur_collect = Vector{Vector{Float64}}(); all_r_ext_collect = Vector{Vector{Float64}}(); 
# all_SA_collect = Vector{Vector{Float64}}(); all_May_collect = Vector{Vector{Float64}}(); 
# all_SAN_collect = Vector{Vector{Float64}}(); all_MayN_collect = Vector{Vector{Float64}}(); 
# all_sym_collect = Vector{Vector{Float64}}()
@time for j in 1: 50
    all_leading_H = ComplexF64[]; all_rich_H = Float64[]; all_r_sur_H = Float64[]; all_r_ext_H = Float64[];
    # all_SA_H = Float64[]; all_May_H = Float64[];
    # all_SAN_H = Float64[]; all_MayN_H = Float64[];
    # all_sym_H = Float64[]
    for i in 1:length(path)
        @load path[i] all_sur all_leading all_r # all_ℵ all_C
        S = length(all_sur[j])
        sur = Int64.(all_sur[j])
        all_r = all_r[j]
        r_sur = all_r[sur]
        r_ext = all_r[setdiff(1:N, sur)]
        # A = reshape(all_ℵ[j], N, N)
        # LV_Jac = [A[i, x]*all_C[j][i] for i in 1:N, x in 1:N]
        # LV_Jac[diagind(LV_Jac)] .= [all_r[j][i] + A[i, i]*all_C[j][i] + sum(A[i, x]*all_C[j][x] for x in 1:N) for i in 1:N]
        # J_ij = [LV_Jac[i, j] for i in 1:N for j in 1:N if i>j]
        # J_ji = [LV_Jac[j, i] for i in 1:N for j in 1:N if i>j]
        # τ = mean(J_ij .* J_ji)*N
        # d = mean(abs.(diag(LV_Jac)))
        # X = [LV_Jac[i,j] for i in 1:N for j in 1:N if i!= j]
        # σ = std(X) 
        # E = mean(X)
        # E2 = mean(J_ij .* J_ji)
        # ρ = (E2-E^2)/(σ^2)
        # θ = d/σ
        # SA = sqrt(S) - θ/(1+E^2/σ^2) # < 0
        # SA_N = sqrt(N) - θ/(1+E^2/σ^2) # < 0
        # May = sqrt(S) * σ - d # < 0
        # May_N = sqrt(N) * σ - d # < 0
        append!(all_leading_H, all_leading[j]); append!(all_rich_H, S);
        append!(all_r_sur_H, r_sur);append!(all_r_ext_H, r_ext);
        # append!(all_SA_H, SA); append!(all_May_H, May);
        # append!(all_SAN_H, SA_N); append!(all_MayN_H, May_N);
        # append!(all_sym_H, τ)
        next!(progress)
    end 
    push!(all_leading_collect, all_leading_H); push!(all_rich_collect, all_rich_H);
    push!(all_r_sur_collect, all_r_sur_H); push!(all_r_ext_collect, all_r_ext_H); 
    # push!(all_SA_collect, all_SA_H); push!(all_May_collect, all_May_H);
    # push!(all_SAN_collect, all_SAN_H); push!(all_MayN_collect, all_MayN_H);
    # push!(all_sym_collect, all_sym_H)
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

sym = [mean(all_sym_collect[v]) for v in 2: 50]
sym[abs.(sym) .> 1.0e3] .= NaN
plot(sym)

SA_sta = [mean(all_SA_collect[v]) for v in 1:50]
May_sta = [mean(all_May_collect[v]) for v in 1:50]

SA_sta_N = [mean(all_SAN_collect[v]) for v in 1:50]
May_sta_N = [mean(all_MayN_collect[v]) for v in 1:50]
May_sta[abs.(May_sta).> 1000] .= NaN
May_sta_N[abs.(May_sta_N).> 1000] .= NaN

# sum(abs.(imag.((vcat(all_leading_collect...)))) .> 0)
plot(May_sta_N, sta_var)
plot(May_sta, sta_var)
plot(SA_sta_N, sta_var)
plot(SA_sta, sta_var)

sta_var = [sum(real.(all_leading_collect[t]) .< 0)/length(path) for t in 1:50]
rich_var = [mean(all_rich_collect[t]) for t in 1: 50]
rich_var_err= [std(all_rich_collect[t])/sqrt(length(all_rich_collect[t])) for t in 1: 50]

using ColorSchemes
cs = ColorScheme(range(colorant"#F8BA17",colorant"#601210",  length = 50))
f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Richness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
sc = scatter!(ax, rich_var, sta_var, color = 1:50, colormap = cs, markersize = 25, alpha = 0.8)
for (x, y, e) in zip(rich_var, sta_var, rich_var_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = (:black, 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(sta_var)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = (:blue, 0.3), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = (:blue, 0.3), linewidth = 1)
end
subgrid = GridLayout(f[1,2], tellheight = false)
# Colorbar(f[1,2], colorrange = [0.0, 1.5], colormap = cs, label = L"σ_{log(B)}", labelsize = 50)
Colorbar(f[1,2], colorrange = [0.0, maximum(var_B)], colormap = cs, label = L"σ_{B}", labelsize = 50)
# Label(f[1,1, TopLeft()], "(b)")
f
# save("../results/var_um_rich_sta.pdf", f) 
save("../results/var_um_rich_sta_withoutlog.pdf", f) 

r_sur_meam = [mean(all_r_sur_collect[t]) for t in 1: 50]
r_sur_meam_err= [std(all_r_sur_collect[t])/sqrt(length(all_r_sur_collect[t])) for t in 1: 50]
r_ext_meam = [mean(all_r_ext_collect[t]) for t in 1: 50]
r_ext_meam_err= [std(all_r_ext_collect[t])/sqrt(length(all_r_ext_collect[t])) for t in 1: 50]
# @load "../data/20240914/v_new/v_$(index).jld2" all_sur all_ℵ all_r all_leading all_diag_dom all_C all_R

##########################################
density(abs.(all_sur_α_collect[10]))
density(abs.(all_α_collect[10]))
all_f = []
for i in 1:10
    f = Figure(fontsize = 35, size = (1200, 900));
    ax = Axis(f[1,1], xlabel = "|α|", ylabel = "density", xlabelsize = 50, ylabelsize = 50)
    density!(ax, abs.(all_α_collect[5*i]), color = ("#376298", 0.4), label = "α")
    density!(ax, abs.(all_sur_α_collect[5*i]), color = ("#9A2B1A", 0.4), label = "survivors")
    axislegend(position = :rt)
    # push!(all_f, f)
    display(f)
end 
