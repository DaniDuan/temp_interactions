include("./sim_frame.jl")
using ProgressMeter, Glob#, ColorSchemes, RCall
num_temps = 31
N=100; M=50
Temp_rich = range(0, num_temps-1, length = num_temps)

############## collecting results ##############
path_1 = glob("Eff_iters0_nol_rho1_*", "../data/20250917/p05_nol/")
progress = Progress(length(path_1)*num_temps; desc="Progress running:")
all_rich_collect_1 = Vector{Vector{Float64}}(); all_ii_collect_1 = Vector{Vector{Float64}}(); all_ij_collect_1 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_rich = Float64[]; α_ij =  Float64[]; α_ii = Float64[];
    for i in 1:length(path_1)
        @load path_1[i]  rich all_ℵii all_ℵij #all_sur rich
        append!(all_rich, rich[j]); append!(α_ii, all_ℵii[j]); append!(α_ij, all_ℵij[j]); 
        next!(progress)
    end 
    push!(all_rich_collect_1, all_rich); push!(all_ii_collect_1, α_ii); push!(all_ij_collect_1, α_ij);
end 
# R"library(beepr); beep(sound = 4, expr = NULL)"

path_2 = glob("Eff_iters0_nol_rho2_*", "../data/20250903/p05_nol/")
progress = Progress(length(path_2)*num_temps; desc="Progress running:")
all_rich_collect_2 = Vector{Vector{Float64}}(); all_ii_collect_2 = Vector{Vector{Float64}}(); all_ij_collect_2 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_rich = Float64[]; α_ij =  Float64[]; α_ii = Float64[];
    for i in 1:length(path_2)
        @load path_2[i]  rich all_ℵii all_ℵij #all_sur rich
        append!(all_rich, rich[j]); append!(α_ii, all_ℵii[j]); append!(α_ij, all_ℵij[j]); 
        next!(progress)
    end 
    push!(all_rich_collect_2, all_rich); push!(all_ii_collect_2, α_ii); push!(all_ij_collect_2, α_ij);
end 
# R"library(beepr); beep(sound = 4, expr = NULL)"

path_3 = glob("Eff_iters0_nol_rho3_*", "../data/20250903/p05_nol/")
progress = Progress(length(path_3)*num_temps; desc="Progress running:")
all_rich_collect_3 = Vector{Vector{Float64}}(); all_ii_collect_3 = Vector{Vector{Float64}}(); all_ij_collect_3 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_rich = Float64[]; α_ij =  Float64[]; α_ii = Float64[];
    for i in 1:length(path_3)
        @load path_3[i]  rich all_ℵii all_ℵij #all_sur rich
        append!(all_rich, rich[j]); append!(α_ii, all_ℵii[j]); append!(α_ij, all_ℵij[j]); 
        next!(progress)
    end 
    push!(all_rich_collect_3, all_rich); push!(all_ii_collect_3, α_ii); push!(all_ij_collect_3, α_ij);
end 
# R"library(beepr); beep(sound = 4, expr = NULL)"

path_4 = glob("Eff_iters0_nol_rho4_*", "../data/20250911/p05_nol/")
progress = Progress(length(path_4)*num_temps; desc="Progress running:")
all_rich_collect_4 = Vector{Vector{Float64}}(); all_ii_collect_4 = Vector{Vector{Float64}}(); all_ij_collect_4 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    all_rich = Float64[]; α_ij =  Float64[]; α_ii = Float64[];
    for i in 1:length(path_4)
        @load path_4[i]  rich all_ℵii all_ℵij #all_sur rich
        append!(all_rich, rich[j]); append!(α_ii, all_ℵii[j]); append!(α_ij, all_ℵij[j]); 
        next!(progress)
    end 
    push!(all_rich_collect_4, all_rich); push!(all_ii_collect_4, α_ii); push!(all_ij_collect_4, α_ij);
end 
# R"library(beepr); beep(sound = 4, expr = NULL)"

# @save "../data/20250911/p05_nol/summary_temp_inter_L0_rho.jld2" all_rich_collect_1 all_ii_collect_1 all_ij_collect_1 all_rich_collect_2 all_ii_collect_2 all_ij_collect_2 all_rich_collect_3 all_ii_collect_3 all_ij_collect_3 all_rich_collect_4 all_ii_collect_4 all_ij_collect_4
# @load "../data/20250903/p05_nol/summary_temp_inter_L0_rho.jld2" all_rich_collect_1 all_ii_collect_1 all_ij_collect_1 all_rich_collect_2 all_ii_collect_2 all_ij_collect_2 all_rich_collect_3 all_ii_collect_3 all_ij_collect_3 all_rich_collect_4 all_ii_collect_4 all_ij_collect_4
# @load "../data/summary_temp_inter_L0_rho.jld2" all_rich_collect_1 all_ii_collect_1 all_ij_collect_1 all_rich_collect_2 all_ii_collect_2 all_ij_collect_2 all_rich_collect_3 all_ii_collect_3 all_ij_collect_3 all_rich_collect_4 all_ii_collect_4 all_ij_collect_4



#################################################################################################################
path_1 = glob("Eff_iters0_nol_rho1_*", "../data/20250917/p05_nol/")
progress = Progress(length(path_1)*num_temps; desc="Progress running:")
all_leading_collect_1 = Vector{Vector{ComplexF64}}(); dℵij_collect_1 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    leading = ComplexF64[]; dℵij_H = Float64[];
    for i in 1:length(path_1)
        @load path_1[i] all_leading all_ℵii all_ℵij 
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(leading, all_leading[j]); append!(dℵij_H, dℵij); 
        next!(progress)
    end 
    push!(all_leading_collect_1, leading); push!(dℵij_collect_1, dℵij_H); 
end 
sta_1 = [sum(real.(all_leading_collect_1[t]) .< 0)/length(path_1) for t in 1:num_temps]

path_2 = glob("Eff_iters0_nol_rho2_*", "../data/20250917/p05_nol/")
progress = Progress(length(path_2)*num_temps; desc="Progress running:")
all_leading_collect_2 = Vector{Vector{ComplexF64}}(); dℵij_collect_2 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    leading = ComplexF64[]; dℵij_H = Float64[];
    for i in 1:length(path_2)
        @load path_2[i] all_leading all_ℵii all_ℵij 
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(leading, all_leading[j]); append!(dℵij_H, dℵij); 
        next!(progress)
    end 
    push!(all_leading_collect_2, leading); push!(dℵij_collect_2, dℵij_H); 
end 
sta_2 = [sum(real.(all_leading_collect_2[t]) .< 0)/length(path_2) for t in 1:num_temps]

path_3 = glob("Eff_iters0_nol_rho3_*", "../data/20250917/p05_nol/")
progress = Progress(length(path_3)*num_temps; desc="Progress running:")
all_leading_collect_3 = Vector{Vector{ComplexF64}}(); dℵij_collect_3 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    leading = ComplexF64[]; dℵij_H = Float64[];
    for i in 1:length(path_3)
        @load path_3[i] all_leading all_ℵii all_ℵij 
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(leading, all_leading[j]); append!(dℵij_H, dℵij); 
        next!(progress)
    end 
    push!(all_leading_collect_3, leading); push!(dℵij_collect_3, dℵij_H); 
end 
sta_3 = [sum(real.(all_leading_collect_3[t]) .< 0)/length(path_3) for t in 1:num_temps]

path_4 = glob("Eff_iters0_nol_rho4_*", "../data/20250917/p05_nol/")
progress = Progress(length(path_4)*num_temps; desc="Progress running:")
all_leading_collect_4 = Vector{Vector{ComplexF64}}(); dℵij_collect_4 = Vector{Vector{Float64}}();
idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i,j] for i in 1:N for j in 1:N if i != j]
@time for j in 1: num_temps
    leading = ComplexF64[]; dℵij_H = Float64[];
    for i in 1:length(path_4)
        @load path_4[i] all_leading all_ℵii all_ℵij 
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]
        dℵij = [A[j, i]/A[j, j] for i in 1:N for j in 1:N if j != i]
        append!(leading, all_leading[j]); append!(dℵij_H, dℵij); 
        next!(progress)
    end 
    push!(all_leading_collect_4, leading); push!(dℵij_collect_4, dℵij_H); 
end 

@save "../data/20250917/p05_nol/summary_stability_L0_rho.jld2" all_leading_collect_1 all_leading_collect_2 all_leading_collect_3 all_leading_collect_4 dℵij_collect_1 dℵij_collect_2 dℵij_collect_3 dℵij_collect_4
