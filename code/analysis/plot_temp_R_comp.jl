include("./sim_frame.jl")
using ProgressMeter, RCall, Glob, ColorSchemes
num_temps = 31
N=100; M=50
Temp_rich = range(0, num_temps-1, length = num_temps)
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]

############## collecting results ##############
path_0 = glob("Eff_iters*", "../data/Eff_p0_Umatrix/")
progress = Progress(length(path_0)*num_temps; desc="Progress running:")
num_temps = 31
dℵij_collect_0 = Vector{Vector{Float64}}() ; all_R_collect_0 = Vector{Vector{Float64}}();
@time for j in 1: num_temps
    dℵij_H = Float64[]; all_R_H = Float64[];
    for i in 1:length(path_0)
        @load path_0[i]  all_ℵii all_ℵij all_R
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(dℵij_H, dℵij); append!(all_R_H, all_R[j]); 
        next!(progress)
    end 
    push!(dℵij_collect_0, dℵij_H); push!(all_R_collect_0, all_R_H);
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

path_1 = glob("Eff_iters*", "../data/Eff_p-1_Umatrix/")
progress = Progress(length(path_1)*num_temps; desc="Progress running:")
num_temps = 31
dℵij_collect_1 = Vector{Vector{Float64}}() ; all_R_collect_1 = Vector{Vector{Float64}}();
@time for j in 1: num_temps
    dℵij_H = Float64[]; all_R_H = Float64[];
    for i in 1:length(path_1)
        @load path_1[i]  all_ℵii all_ℵij all_R
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(dℵij_H, dℵij); append!(all_R_H, all_R[j]); 
        next!(progress)
    end 
    push!(dℵij_collect_1, dℵij_H); push!(all_R_collect_1, all_R_H);
end 
R"library(beepr); beep(sound = 4, expr = NULL)"

############## analysing ##############
αij_d_0 = [mean(dℵij_collect_0[t]) for t in 1:num_temps]
αij_d_err_0 = [std(dℵij_collect_0[t])/sqrt(length(dℵij_collect_0[t])) for t in 1: num_temps]
meanR_0 = [mean(all_R_collect_0[t]) for t in 1:num_temps]
R_err_0 = [std(all_R_collect_0[t])/sqrt(length(all_R_collect_0[t])) for t in 1: num_temps]

αij_d_1 = [mean(dℵij_collect_1[t]) for t in 1:num_temps]
αij_d_err_1 = [std(dℵij_collect_1[t])/sqrt(length(dℵij_collect_1[t])) for t in 1: num_temps]
meanR_1 = [mean(all_R_collect_1[t]) for t in 1:num_temps]
R_err_1 = [std(all_R_collect_1[t])/sqrt(length(all_R_collect_1[t])) for t in 1: num_temps]

############## plotting ##############
f = Figure(fontsize = 35, size = (2400, 900));
# Label(f[0,1], "Minimal Trade-off", fontsize = 50)
ax1 = Axis(f[1,1], title = "Minimal Trade-off", xlabel = "Temperature (°C)", ylabel = "Effective Resource Competition\n(αji/αii)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "Resource Abundance", yaxisposition = :right, yticklabelalign = (:left, :center), xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
# hidedecorations!(ax3, grid = false, ticks = true, ticklabels = true)
lines!(ax1, Temp_rich, αij_d_0, color = ("#376298",0.8), linewidth = 5, label = "")
band!(ax1, Temp_rich, αij_d_0 .- αij_d_err_0 , αij_d_0.+ αij_d_err_0, color = ("#376298", 0.3))
lines!(ax2, Temp_rich, meanR_0, color =( "#F8BA17", 0.9), linewidth = 5, label = "")
band!(ax2, Temp_rich, meanR_0 .- R_err_0 , meanR_0 .+ R_err_0 , color = ("#F8BA17", 0.5))
linkxaxes!(ax1,ax2)
lines!(ax1, [0, 30], [1, 1], linestyle = :dash, color = ("#4F363E", 0.9), linewidth = 2)
text!(ax1, 0, 1.05, text = "↑ αⱼᵢ > αᵢᵢ", align = (:left, :center),fontsize = 30)
text!(ax1, 0, 0.95, text = "↓ αⱼᵢ < αᵢᵢ", align = (:left, :center),fontsize = 30)
l1 = [LineElement(color = ("#376298", 0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#F8BA17", 0.9), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, [ "αji/αii", "Resource"], halign = :center, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,1, TopLeft()], "(a)")

ax3 = Axis(f[1,2], title = "Maximal Trade-off", xlabel = "Temperature (°C)", ylabel = "Effective Resource Competition\n(αji/αii)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax4 = Axis(f[1,2], ylabel = "Resource Abundance", yaxisposition = :right, yticklabelalign = (:left, :center), xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
# hidedecorations!(ax3, grid = false, ticks = true, ticklabels = true)
lines!(ax3, Temp_rich, αij_d_1, color = ("#376298",0.8), linewidth = 5, label = "")
band!(ax3, Temp_rich, αij_d_1 .- αij_d_err_1 , αij_d_1.+ αij_d_err_1, color = ("#376298", 0.3))
lines!(ax4, Temp_rich, meanR_1, color =( "#F8BA17", 0.9), linewidth = 5, label = "")
band!(ax4, Temp_rich, meanR_1 .- R_err_1 , meanR_1 .+ R_err_1 , color = ("#F8BA17", 0.5))
linkxaxes!(ax1,ax2)
lines!(ax3, [0, 30], [1, 1], linestyle = :dash, color = ("#4F363E", 0.9), linewidth = 2)
text!(ax3, 0, 1.05, text = "↑ αⱼᵢ > αᵢᵢ", align = (:left, :center),fontsize = 30)
text!(ax3, 0, 0.95, text = "↓ αⱼᵢ < αᵢᵢ", align = (:left, :center),fontsize = 30)
Legend(f[1,2], [l1, l2], tellheight = false, tellwidth = false, [ "αji/αii", "Resource"], halign = :center, valign = :top, framevisible = false) # "ƒc-ƒo"
Label(f[1,2, TopLeft()], "(b)")
f

save("../results/αR.pdf", f) 

# ########### 
# path_0 = glob("Eff_iters*", "../data/Eff_p0_Umatrix/")
# progress = Progress(length(path_0)*num_temps; desc="Progress running:")
# num_temps = 31
# all_sur_α_collect_0 = Vector{Vector{Float64}}(); all_α_collect_0 = Vector{Vector{Float64}}();
# σ_u_collect_0 = Vector{Vector{Float64}}(); mean_α_collect_0 = Vector{Vector{Float64}}(); mean_sur_α_collect_0 = Vector{Vector{Float64}}(); 
# err_α_collect_0 = Vector{Vector{Float64}}(); err_sur_α_collect_0 = Vector{Vector{Float64}}(); 
# # all_αijii_collect_0 = Vector{Vector{Float64}}()
# idx = collect(CartesianIndices(zeros(Float64, N, N)))
# ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
# @time for j in 1: num_temps
#     all_sur_α = Float64[]; all_α = Float64[];
#     σ_u_H = Float64[]; mean_α_H = Float64[]; mean_sur_α_H = Float64[]; 
#     err_α_H = Float64[]; err_sur_α_H = Float64[]; 
#     for i in 1:length(path_0)
#         @load path_0[i]  all_ℵii all_ℵii_sur all_ℵij all_R all_u
#         sur = findall(x -> x in all_ℵii_sur[j], all_ℵii[j])
#         A = zeros(Float64, N, N)
#         A[ind_off] = all_ℵij[j]
#         A[diagind(A)] = all_ℵii[j]
#         α_sur = A[sur, sur]
#         α_sur_off = [α_sur[i,j] for i in 1:length(sur) for j in 1:length(sur) if i != j]
#         dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]

#         σ_u = std(log.(abs.(all_u[j]))); mean_α = mean(log.(abs.(all_ℵij[j]))); mean_sur_α = mean(log.(abs.(α_sur_off)));
#         err_α = std(log.(abs.(all_ℵij[j])))/sqrt(length(all_ℵij[j])); err_sur_α = std(log.(abs.(α_sur_off)))/sqrt(length(α_sur_off)); 
#         append!(all_sur_α, α_sur_off); append!(all_α, all_ℵij[j]);
#         append!(σ_u_H, σ_u); append!(mean_α_H, mean_α); append!(mean_sur_α_H, mean_sur_α); 
#         append!(err_α_H, err_α); append!(err_sur_α_H, err_sur_α); 
#         next!(progress)
#     end 
#     push!(all_sur_α_collect_0, all_sur_α); push!(all_α_collect_0, all_α);
#     push!(σ_u_collect_0, σ_u_H); push!(mean_α_collect_0, mean_α_H); push!(mean_sur_α_collect_0, mean_sur_α_H); 
#     push!(err_α_collect_0, err_α_H); push!(err_sur_α_collect_0, err_sur_α_H); 
# end 
# R"library(beepr); beep(sound = 4, expr = NULL)"

# f = Figure(fontsize = 35, size = (1200, 900));
# ax = Axis(f[1,1], xlabel = "|α|", ylabel = "density", xlabelsize = 50, ylabelsize = 50)
# density!(ax, abs.(all_α_collect[10]), color = ("#376298", 0.4), label = "α")
# density!(ax, abs.(all_sur_α_collect[10]), color = ("#9A2B1A", 0.4), label = "survivors")
# axislegend(position = :rt)
# # push!(all_f, f)
# display(f)

# α_sur_mean = [mean(log.(abs.(all_sur_α_collect_0[t]))) for t in 1: num_temps]
# α_sur_err= [std(log.(abs.(all_sur_α_collect_0[t])))/sqrt(length(all_sur_α_collect_0[t])) for t in 1: num_temps]
# α_mean = [mean(log.(abs.(all_α_collect_0[t]))) for t in 1: num_temps]
# α_err= [std(log.(abs.(all_α_collect_0[t])))/sqrt(length(all_α_collect_0[t])) for t in 1: num_temps]

# f = Figure(fontsize = 30, size = (1200, 900));
# ax1 = Axis(f[1,1], title = "Minimal Trade-off", xlabel = "Temperature", ylabel = "log(|α|)", xlabelsize = 35, ylabelsize = 35, ygridvisible = false, xgridvisible = false)
# # xlims!(nothing, 7)
# lines!(ax1, Temp_rich,α_sur_mean, color = ("#601210", 0.7), linewidth = 5, label = "survivors")
# band!(ax1, Temp_rich, α_sur_mean .- α_sur_err, α_sur_mean .+ α_sur_err, color = ("#601210", 0.2))
# lines!(ax1, Temp_rich,α_mean, color = ("#5676A5", 0.7), linewidth = 5, label = "α")
# band!(ax1, Temp_rich, α_mean .- α_err, α_mean .+ α_err, color = ("#5676A5", 0.2))
# axislegend(position = :rt)
# # Label(f[1,1, TopLeft()], "(c)")
# f
# save("../results/α_sur0.pdf", f) 

# σ_u = vcat(σ_u_collect...)
# all_α = vcat(mean_α_collect...)
# all_sur_α = vcat(mean_sur_α_collect...)
# all_α_err = vcat(err_α_collect...)
# all_sur_α_err = vcat(err_sur_α_collect...)

# f = Figure(fontsize = 35, size = (1200, 900));
# ax = Axis(f[1,1], title = "Minimal Trade-off", xlabel = "σ(log(u))", ylabel = "log(|α|)", xlabelsize = 50, ylabelsize = 50)
# scatter!(ax, σ_u, all_sur_α, color = "#601210", markersize = 10, alpha = 0.5, label = "survivors")
# for (x, y, e) in zip(σ_u, all_sur_α, all_sur_α_err)
#     # Vertical line
#     lines!(ax, [x, x], [y - e, y + e], color = ("#601210", 0.3), linewidth = 1)
#     # lines!(ax, [x - e, x + e], [y, y], color = ("#601210", 0.4), linewidth = 1)
#     # Horizontal caps
#     cap_length = 0.0005 * mean(all_sur_α)  # Length of horizontal caps
#     lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#601210", 0.3), linewidth = 1)
#     lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#601210", 0.3), linewidth = 1)
#     # lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = ("#601210", 0.4), linewidth = 1)
#     # lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = ("#601210", 0.4), linewidth = 1)
# end
# scatter!(ax, σ_u, all_α, color = "#5676A5", markersize = 10, alpha = 0.5, label = "α")
# for (x, y, e) in zip(σ_u, all_α, all_α_err)
#     # Vertical line
#     lines!(ax, [x, x], [y - e, y + e], color = ("#5676A5", 0.3), linewidth = 1)
#     # lines!(ax, [x - e, x + e], [y, y], color = ("#5676A5", 0.4), linewidth = 1)
#     # Horizontal caps
#     cap_length = 0.0005 * mean(all_α)  # Length of horizontal caps
#     lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#5676A5", 0.3), linewidth = 1)
#     lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#5676A5", 0.3), linewidth = 1)
#     # lines!(ax, [x - e, x - e], [y - cap_length, y + cap_length], color = ("#5676A5", 0.4), linewidth = 1)
#     # lines!(ax, [x + e, x + e], [y - cap_length, y + cap_length], color = ("#5676A5", 0.4), linewidth = 1)
# end
# axislegend(position = :rb)
# # Label(f[1,1, TopLeft()], "(a)")
# f
# save("../results/α_sur0_scatter.pdf", f) 


