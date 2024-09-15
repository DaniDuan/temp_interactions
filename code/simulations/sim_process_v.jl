include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob

N=100
M=50
### Temp params 
# ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
Ci = fill(0.1, N)
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator, N, M) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

path = glob("v_*", "../data/v_new/")
progress = Progress(length(path)*50; desc="Progress running:")
num_temps = 31
all_leading_collect = Vector{Vector{ComplexF64}}(); all_rich_collect = Vector{Vector{Float64}}()
@time for j in 1: 50
    all_leading_H = ComplexF64[]; all_rich_H = Float64[]
    for i in 1:length(path)
        @load path[i] all_sur all_ℵ all_r all_leading all_diag_dom
        append!(all_leading_H, all_leading[j]);
        append!(all_rich_H, length(all_sur[j]))
        next!(progress)
    end 
    push!(all_leading_collect, all_leading_H); push!(all_rich_collect, all_rich_H)
end 

# sum(abs.(imag.((vcat(all_leading_collect...)))) .> 0)

sta_var = [sum(real.(all_leading_collect[t]) .< 0)/length(path) for t in 1:50]
rich_var = [mean(all_rich_collect[t]) for t in 1: 50]
rich_var_err= [std(all_rich_collect[t])/sqrt(length(all_rich_collect[t])) for t in 1: 50]

using ColorSchemes
cs = ColorScheme(range(colorant"#376298",colorant"#9A2B1A",  length = 50))
f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Richness", ylabel = "p(Stability)", xlabelsize = 50, ylabelsize = 50)
sc = scatter!(ax, rich_var, sta_var, color = 1:50, colormap = cs, markersize = 15, alpha = 0.8)
for (x, y, e) in zip(rich_var, sta_var, rich_var_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color = (:black, 0.4), linewidth = 1)
    # Horizontal caps
    cap_length = 0.001 * mean(sta_var)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = (:black, 0.4), linewidth = 1)
    lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = (:black, 0.4), linewidth = 1)
end
subgrid = GridLayout(f[1,2], tellheight = false)
Colorbar(f[1,2], colorrange = [0.0, 1.5], colormap = cs, label = "σ")
f
save("../results/var_um_rich_sta.pdf", f) 
