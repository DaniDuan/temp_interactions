include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob
using ColorSchemes

###### richness 
rich_L = zeros(7, num_temps); rich_Lerr = zeros(7, num_temps)
progress = Progress(length(path)*num_temps*7; desc="Progress running:")
for leakage in 1:7
    path = glob("*L$(leakage)*", "../data/L_re/")
    num_temps = 31
    rich_collect = Vector{Vector{Float64}}()#; all_ij_collect = Vector{Vector{Float64}}()
    @time for j in 1: num_temps
        rich_H = Float64[]#; all_ij_H = Float64[]
        for i in 1:length(path)
            @load path[i]  rich# all_ℵij #all_ℵii_sur
            push!(rich_H, rich[j])#; append!(all_ij_H, all_ℵij[j])
            next!(progress)
        end 
        push!(rich_collect, rich_H)#; push!(all_ij_collect, all_ij_H)
    end 
    rich_L[leakage,:] = [mean(rich_collect[t]) for t in 1:num_temps]
    rich_Lerr[leakage,:] = [std(rich_collect[t])/sqrt(length(rich_collect[t])) for t in 1:num_temps]
end 

Temp_rich = range(0, num_temps-1, length = num_temps)

num_temps = 31
# cs = ColorScheme(range(colorant"#FA8328",colorant"#015845", length = 7))
cscheme = ColorScheme(range(colorant"#376298",colorant"#F5D44B", length = 4))
cscheme1 = ColorScheme(range(colorant"#F5D44B",colorant"#9A2B1A", length = 4))
cs = vcat(cscheme[1:4], cscheme1[2:4])

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Richness", xlabelsize = 50, ylabelsize = 50)
for leakage in 1:7
    lines!(ax, Temp_rich,rich_L[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
    band!(ax, Temp_rich, rich_L[leakage,:] .- rich_Lerr[leakage,:], rich_L[leakage,:] .+ rich_Lerr[leakage,:], color = (cs[leakage], 0.3))
end 
Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
# Label(f[1,1, TopLeft()], "ρ = 0")
Label(f[1,1, TopLeft()], "realistic ρ")
f
save("../results/rich_L_re.pdf", f) 


########### interactions 
mean_ii = zeros(7, num_temps); mean_iierr = zeros(7, num_temps);
mean_ij = zeros(7, num_temps); mean_ijerr = zeros(7, num_temps);
mean_d = zeros(7, num_temps); mean_derr = zeros(7, num_temps);
# progress = Progress(length(path)*num_temps*7; desc="Progress running:")
for leakage in 1:7
    path = glob("*L$(leakage)*", "../data/L0/")
    num_temps = 31
    all_ii_collect = Vector{Vector{Float64}}(); all_ij_collect = Vector{Vector{Float64}}(); all_ℵij_d_collect = Vector{Vector{Float64}}()
    @time for j in 1: num_temps
        all_ii_H = Float64[]; all_ij_H = Float64[]; all_ℵij_d_H = Float64[]
        for i in 1:length(path)
            @load path[i] all_ℵii all_ℵij all_ℵij_d
            append!(all_ii_H, all_ℵii[j]); append!(all_ij_H, all_ℵij[j]); append!(all_ℵij_d_H, all_ℵij_d[j])
            # next!(progress)
        end 
        push!(all_ii_collect, all_ii_H); push!(all_ij_collect, all_ij_H); push!(all_ℵij_d_collect, all_ℵij_d_H)
    end 
    mean_ii[leakage,:] = [mean(all_ii_collect[t]) for t in 1:num_temps]
    mean_iierr[leakage,:] = [std(all_ii_collect[t])/sqrt(length(all_ii_collect[t])) for t in 1:num_temps]
    mean_ij[leakage,:] = [mean(all_ij_collect[t]) for t in 1:num_temps]
    mean_ijerr[leakage,:] = [std(all_ij_collect[t])/sqrt(length(all_ij_collect[t])) for t in 1:num_temps]
    mean_d[leakage,:] = [mean(all_ℵij_d_collect[t]) for t in 1:num_temps]
    mean_derr[leakage,:] = [std(all_ℵij_d_collect[t])/sqrt(length(all_ℵij_d_collect[t])) for t in 1:num_temps]
end 

num_temps = 31
Temp_rich = range(0, num_temps-1, length = num_temps)

# cs = ColorScheme(range(colorant"#FA8328",colorant"#015845", length = 7))
cscheme = ColorScheme(range(colorant"#376298",colorant"#F5D44B", length = 4))
cscheme1 = ColorScheme(range(colorant"#F5D44B",colorant"#9A2B1A", length = 4))
cs = vcat(cscheme[1:4], cscheme1[2:4])

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], title = "Minimal Trade-off", xlabel = "Temperature (°C)", ylabel = L"α_{ii}", xlabelsize = 50, ylabelsize = 50)
for leakage in 1:7
    lines!(ax, Temp_rich,mean_ii[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
    band!(ax, Temp_rich, mean_ii[leakage,:] .- mean_iierr[leakage,:], mean_ii[leakage,:] .+ mean_iierr[leakage,:], color = (cs[leakage], 0.3))
end 
Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
# Label(f[1,1, TopLeft()], (a))
# Label(f[1,1, TopLeft()], "realistic ρ")
f
save("../results/ii_L0.pdf", f) 

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], title = "Minimal Trade-off", xlabel = "Temperature (°C)", ylabel = L"α_{i≠j}", xlabelsize = 50, ylabelsize = 50)
for leakage in 1:7
    lines!(ax, Temp_rich,mean_ij[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
    band!(ax, Temp_rich, mean_ij[leakage,:] .- mean_ijerr[leakage,:], mean_ij[leakage,:] .+ mean_ijerr[leakage,:], color = (cs[leakage], 0.3))
end 
Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
# Label(f[1,1, TopLeft()], "ρ = -1")
# Label(f[1,1, TopLeft()], "realistic ρ")
f
save("../results/ij_L0.pdf", f) 

# f = Figure(fontsize = 35, size = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αᵢᵢ", xlabelsize = 50, ylabelsize = 50)
# for leakage in 1:7
#     lines!(ax, Temp_rich,mean_d[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
#     band!(ax, Temp_rich, mean_d[leakage,:] .- mean_derr[leakage,:], mean_d[leakage,:] .+ mean_derr[leakage,:], color = (cs[leakage], 0.3))
# end 
# Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
# Label(f[1,1, TopLeft()], "ρ = 0")
# # Label(f[1,1, TopLeft()], "realistic ρ")
# f
# save("../results/ij_ii_L0.pdf", f) 

########### 
ρ_ijji = Float64[]
progress = Progress(length(path)*num_temps*7; desc="Progress running:")
for leakage in 1:7
    path = glob("*L$(leakage)*", "../data/L0/")
    num_temps = 31
    all_uℵij_collect = Vector{Vector{Float64}}(); all_lℵij_collect = Vector{Vector{Float64}}()
    @time for j in 1: num_temps
        all_uℵij_H = Float64[]; all_lℵij_H = Float64[]
        for i in 1:length(path)
            @load path[i] all_uℵij all_lℵij
            append!(all_uℵij_H, all_uℵij[j]); append!(all_lℵij_H, all_lℵij[j])
            next!(progress)
        end 
        push!(all_uℵij_collect, all_uℵij_H); push!(all_lℵij_collect, all_lℵij_H)
    end 
    push!(ρ_ijji, cor(vcat(all_uℵij_collect...), vcat(all_lℵij_collect...)))
end 
all_ρα =  DataFrame(Leakage = range(0.1,0.7,7), ρ_α = ρ_ijji)
num_temps = 31
# cs = ColorScheme(range(colorant"#FA8328",colorant"#015845", length = 7))
cscheme = ColorScheme(range(colorant"#376298",colorant"#F5D44B", length = 4))
cscheme1 = ColorScheme(range(colorant"#F5D44B",colorant"#9A2B1A", length = 4))
cs = vcat(cscheme[1:4], cscheme1[2:4])

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Richness", xlabelsize = 50, ylabelsize = 50)
for leakage in 1:7
    lines!(ax, Temp_rich,rich_L[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
    band!(ax, Temp_rich, rich_L[leakage,:] .- rich_Lerr[leakage,:], rich_L[leakage,:] .+ rich_Lerr[leakage,:], color = (cs[leakage], 0.3))
end 
Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
# Label(f[1,1, TopLeft()], "ρ = 0")
Label(f[1,1, TopLeft()], "realistic ρ")
f


################### niche_diff ################################

path = glob("*n3*", "../data/N-1/")
progress = Progress(length(path)*num_temps; desc="Progress running:")
num_temps = 31
rich_collect_1 = Vector{Vector{Float64}}(); all_ℵii_collect_1 = Vector{Vector{Float64}}(); all_ℵij_collect_1 = Vector{Vector{Float64}}()
@time for j in 1: num_temps
    rich_H = Float64[]; all_ℵii_H = Float64[]; all_ℵij_H = Float64[]
    for i in 1:length(path)
        @load path[i]  rich all_ℵii all_ℵij
        push!(rich_H, rich[j]); append!(all_ℵii_H, all_ℵii[j]); append!(all_ℵij_H, all_ℵij[j])
        next!(progress)
    end 
    push!(rich_collect_1, rich_H); push!(all_ℵii_collect_1, all_ℵii_H); push!(all_ℵij_collect_1, all_ℵij_H)
end 
# rich_N[niche-1,:] = [mean(rich_collect[t]) for t in 1:num_temps]
# rich_Nerr[niche-1,:] = [std(rich_collect[t])/sqrt(length(rich_collect[t])) for t in 1:num_temps]
mean_rich = [mean(rich_collect[t]) for t in 1:num_temps]
rich_err = [std(rich_collect[t])/sqrt(length(rich_collect[t])) for t in 1:num_temps]

num_temps = 31
Temp_rich = range(0, num_temps-1, length = num_temps)

f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], title = "Minimum Overlap", xlabel = "Temperature (°C)", ylabel = "Richness", xlabelsize = 50, ylabelsize = 50)
lines!(ax1, Temp_rich, mean_rich, color = ("#FA8328", 0.8), linewidth = 5)
band!(ax1, Temp_rich, mean_rich .- rich_err, mean_rich .+ rich_err, color = ("#FA8328 ", 0.3))
f
save("../results/rich_N0.pdf", f) 



N = 100; M = 50
niche_diff = ones(M, N)
rand_pos = randperm(M)
rand_pos_2 = randperm(M)
for i in 1:N 
    if i <= M
        i_con = rand_pos[i]
    else 
        i_con = rand_pos_2[i - M]
    end 
    niche_diff[i_con, i] = 100.0
end 
u_diff = fu_niche(N, M, niche_diff) # Completely differentiated 

ii_mean_0 =  [mean(all_ℵii_collect_0[t]) for t in 1:num_temps]
ii_err_0 = [std(all_ℵii_collect_0[t])/sqrt(length(all_ℵii_collect_0[t])) for t in 1:num_temps]
ij_mean_0 = [mean(all_ℵij_collect_0[t]) for t in 1:num_temps]
ij_err_0 = [std(all_ℵij_collect_0[t])/sqrt(length(all_ℵij_collect_0[t])) for t in 1:num_temps]
ii_mean_1 =  [mean(all_ℵii_collect_1[t]) for t in 1:num_temps]
ii_err_1 = [std(all_ℵii_collect_1[t])/sqrt(length(all_ℵii_collect_1[t])) for t in 1:num_temps]
ij_mean_1 = [mean(all_ℵij_collect_1[t]) for t in 1:num_temps]
ij_err_1 = [std(all_ℵij_collect_1[t])/sqrt(length(all_ℵij_collect_1[t])) for t in 1:num_temps]


f = Figure(fontsize = 35, size = (2100, 700));
ax1 = Axis(f[1,1][1,1], title = "High niche differentiation", ygridvisible = false, xgridvisible = false, xlabelsize = 35, ylabelsize = 35)
hm1 = heatmap!(ax1, u_diff, colormap = Reverse(:grayC))
ax1.xticks = ([], []); ax1.yticks = ([], [])
subgrid = GridLayout(f[1,1][1,2], tellheight = false)
Label(subgrid[1,1], "High")
Colorbar(subgrid[2,1], hm1, ticksvisible = false, ticklabelsvisible = false)
Label(subgrid[3,1], "Low")
ax2 = Axis(f[1,2], title = "Minimal Trade-off",
    xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax2, Temp_rich, ii_mean_0, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax2, Temp_rich, ii_mean_0 .- ii_err_0, ii_mean_0 .+ ii_err_0, color = ("#FA8328", 0.2))
lines!(ax2, Temp_rich, ij_mean_0, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i ≠ j}")
band!(ax2, Temp_rich,  ij_mean_0 .- ij_err_0, ij_mean_0 .+ ij_err_0, color = ("#015845", 0.2))
axislegend(position = :lb)
ax3 = Axis(f[1,3], title = "Maximal Trade-off",
    xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax3, Temp_rich, ii_mean_1, color = ("#FA8328", 0.8), linewidth = 5, label = L"α_{ii}")
band!(ax3, Temp_rich, ii_mean_1 .- ii_err_1, ii_mean_1 .+ ii_err_1, color = ("#FA8328", 0.2))
lines!(ax3, Temp_rich, ij_mean_1, color = ("#015845", 0.8), linewidth = 5, label = L"α_{i ≠ j}")
band!(ax3, Temp_rich,  ij_mean_1 .- ij_err_1, ij_mean_1 .+ ij_err_1, color = ("#015845", 0.2))
axislegend(position = :lb)
f
save("../results/α_Niche.pdf", f) 





f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], title = "Minimum Overlap",
    xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
lines!(ax1, Temp_rich, ii_mean, color = ("#FA8328", 0.8), linewidth = 5, label = "αii")
band!(ax1, Temp_rich, ii_mean .- ii_err, ii_mean .+ ii_err, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, ij_mean, color = ("#015845", 0.8), linewidth = 5, label = "αij")
band!(ax1, Temp_rich,  ij_mean .- ij_err, ij_mean .+ ij_err, color = ("#015845", 0.2))
axislegend(position = :lb)
# Label(f[1,1, TopLeft()], "(d)")
f
