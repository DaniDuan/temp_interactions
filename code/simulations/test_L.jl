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
progress = Progress(length(path)*num_temps*7; desc="Progress running:")
for leakage in 1:7
    path = glob("*L$(leakage)*", "../data/L0/")
    num_temps = 31
    all_ii_collect = Vector{Vector{Float64}}(); all_ij_collect = Vector{Vector{Float64}}(); all_ℵij_d_collect = Vector{Vector{Float64}}()
    @time for j in 1: num_temps
        all_ii_H = Float64[]; all_ij_H = Float64[]; all_ℵij_d_H = Float64[]
        for i in 1:length(path)
            @load path[i] all_ℵii all_ℵij all_ℵij_d
            append!(all_ii_H, all_ℵii[j]); append!(all_ij_H, all_ℵij[j]); append!(all_ℵij_d_H, all_ℵij_d[j])
            next!(progress)
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

Temp_rich = range(0, num_temps-1, length = num_temps)

num_temps = 31
# cs = ColorScheme(range(colorant"#FA8328",colorant"#015845", length = 7))
cscheme = ColorScheme(range(colorant"#376298",colorant"#F5D44B", length = 4))
cscheme1 = ColorScheme(range(colorant"#F5D44B",colorant"#9A2B1A", length = 4))
cs = vcat(cscheme[1:4], cscheme1[2:4])

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αᵢᵢ", xlabelsize = 50, ylabelsize = 50)
for leakage in 1:7
    lines!(ax, Temp_rich,mean_ii[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
    band!(ax, Temp_rich, mean_ii[leakage,:] .- mean_iierr[leakage,:], mean_ii[leakage,:] .+ mean_iierr[leakage,:], color = (cs[leakage], 0.3))
end 
Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
Label(f[1,1, TopLeft()], "ρ = 0")
# Label(f[1,1, TopLeft()], "realistic ρ")
f
save("../results/ii_L0.pdf", f) 

f = Figure(fontsize = 35, size = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "αᵢⱼ", xlabelsize = 50, ylabelsize = 50)
for leakage in 1:7
    lines!(ax, Temp_rich,mean_ij[leakage,:], color = (cs[leakage], 0.8), linewidth = 5, label = "L = $(leakage * 0.1)")
    band!(ax, Temp_rich, mean_ij[leakage,:] .- mean_ijerr[leakage,:], mean_ij[leakage,:] .+ mean_ijerr[leakage,:], color = (cs[leakage], 0.3))
end 
Colorbar(f[1,2], colorrange = [0.1, 0.7], colormap = cs, label = "Leakage")
Label(f[1,1, TopLeft()], "ρ = 0")
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
