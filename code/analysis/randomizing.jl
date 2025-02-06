include("./sim_frame.jl")
using ProgressMeter, RCall
using Glob
using ColorSchemes
path = glob("Eff_iters*", "../data/Eff_p-1_new/")
# path = glob("Eff_iters*", "../data/L07/p-1/")

N = 100
M = 50
### Temp params 
# ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr = 273.15 + 10;
Ed = 3.5;
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 31

Temp_rich = range(0, num_temps - 1, length=num_temps)

path_1 = glob("Eff_iters*", "../data/Eff_p-1_new/")
progress = Progress(length(path) * num_temps; desc="Progress running:")
all_richness_s0_collect_1 = Vector{Vector{Float64}}();
all_leading_s0_collect_1 = Vector{Vector{ComplexF64}}();
αii_var_s0_collect_1 = Vector{Vector{Float64}}();
αij_var_s0_collect_1 = Vector{Vector{Float64}}();
all_richness_s1_collect_1 = Vector{Vector{Float64}}();
all_leading_s1_collect_1 = Vector{Vector{ComplexF64}}();
αii_var_s1_collect_1 = Vector{Vector{Float64}}();
αij_var_s1_collect_1 = Vector{Vector{Float64}}();
all_richness_s2_collect_1 = Vector{Vector{Float64}}();
all_leading_s2_collect_1 = Vector{Vector{ComplexF64}}();
αii_var_s2_collect_1 = Vector{Vector{Float64}}();
αij_var_s2_collect_1 = Vector{Vector{Float64}}();
all_richness_s3_collect_1 = Vector{Vector{Float64}}();
all_leading_s3_collect_1 = Vector{Vector{ComplexF64}}();
αii_var_s3_collect_1 = Vector{Vector{Float64}}();
αij_var_s3_collect_1 = Vector{Vector{Float64}}();

idx = collect(CartesianIndices(zeros(Float64, N, N)))
ind_off = [idx[i, j] for i in 1:N for j in 1:N if i != j]
@time for j in 1:num_temps
    rich_s0 = Float64[]
    all_leading_s0 = ComplexF64[]
    αii_var_s0 = Float64[]
    αij_var_s0 = Float64[]
    rich_s1 = Float64[]
    all_leading_s1 = ComplexF64[]
    αii_var_s1 = Float64[]
    αij_var_s1 = Float64[]
    rich_s2 = Float64[]
    all_leading_s2 = ComplexF64[]
    αii_var_s2 = Float64[]
    αij_var_s2 = Float64[]
    rich_s3 = Float64[]
    all_leading_s3 = ComplexF64[]
    αii_var_s3 = Float64[]
    αij_var_s3 = Float64[]
    for i in 1:length(path)
        @load path[i] all_ℵii all_ℵij all_r all_C all_ℵii_sur
        sur = findall(x -> x in all_ℵii_sur[j], all_ℵii[j])
        N_s = length(sur)
        C = all_C[j]
        r = all_r[j]
        Ci = fill(0.1, N_s)
        A = zeros(Float64, N, N)
        A[ind_off] = all_ℵij[j]
        A[diagind(A)] = all_ℵii[j]

        # scenario 0: the original
        A_sur = A[sur, sur]
        A_sur_off = [A_sur[i, j] for i in 1:N_s for j in 1:N_s if i != j]
        LV_Jac_s0 = [A_sur[i, j] * C[i] for i in 1:N_s, j in 1:N_s]
        jac_eigen_s0 = eigen(LV_Jac_s0).values
        leading_s0 = jac_eigen_s0[argmax(real.(jac_eigen_s0))]
        append!(rich_s0, N_s)
        append!(all_leading_s0, leading_s0)
        append!(αii_var_s0, var(log.(abs.(diag(A_sur)))))
        append!(αij_var_s0, var(log.(abs.(A_sur_off))))

        # scenario 1: shuffling diagonal
        A_s1 = copy(A_sur)
        for perm in 1:10
            A_s1[diagind(A_s1)] = shuffle(diag(A_sur))
            p_lv_s1 = (ℵ=A_s1, r=r[sur], N=N_s)
            prob_LV_s1 = ODEProblem(LV_dx!, Ci, tspan, p_lv_s1)
            sol_LV_s1 = solve(prob_LV_s1, AutoVern7(Rodas5()), save_everystep=false, callback=cb)
            bm_LV_s1 = sol_LV_s1.u[length(sol_LV_s1.t)][1:N_s]
            N_LV_s1 = length((1:N_s)[bm_LV_s1.>1.0e-7])

            A_s1_off = [A_s1[i, j] for i in 1:N_s for j in 1:N_s if i != j]
            Jac_s1 = [A_s1[i, j] * C[i] for i in 1:N_s, j in 1:N_s]
            jac_eigen_s1 = eigen(Jac_s1).values
            leading_s1 = jac_eigen_s1[argmax(real.(jac_eigen_s1))]
            append!(rich_s1, N_LV_s1)
            append!(all_leading_s1, leading_s1)
            append!(αii_var_s1, var(log.(abs.(diag(A_s1)))))
            append!(αij_var_s1, var(log.(abs.(A_s1_off))))
        end

        # scenario 2: setting diagonal to a constant mean value
        A_s2 = copy(A_sur)
        A_s2[diagind(A_s2)] = fill(-exp(mean(log.(abs.(diag(A_sur))))), N_s)
        p_lv_s2 = (ℵ=A_s2, r=r[sur], N=N_s)
        prob_LV_s2 = ODEProblem(LV_dx!, Ci, tspan, p_lv_s2)
        sol_LV_s2 = solve(prob_LV_s2, AutoVern7(Rodas5()), save_everystep=false, callback=cb)
        bm_LV_s2 = sol_LV_s2.u[length(sol_LV_s2.t)][1:N_s]
        N_LV_s2 = length((1:N_s)[bm_LV_s2.>1.0e-7])

        A_s2_off = [A_s2[i, j] for i in 1:N_s for j in 1:N_s if i != j]
        Jac_s2 = [A_s2[i, j] * C[i] for i in 1:N_s, j in 1:N_s]
        jac_eigen_s2 = eigen(Jac_s2).values
        leading_s2 = jac_eigen_s2[argmax(real.(jac_eigen_s2))]
        append!(rich_s2, N_LV_s2)
        append!(all_leading_s2, leading_s2)
        append!(αii_var_s2, var(log.(abs.(diag(A_s2)))))
        append!(αij_var_s2, var(log.(abs.(A_s2_off))))

        # scenario 3: shuffling off-diagonal
        A_s3 = copy(A_sur)
        idx_s3 = collect(CartesianIndices(zeros(Float64, N_s, N_s)))
        ind_off_s3 = [idx_s3[i, j] for i in 1:N_s for j in 1:N_s if i != j]
        for perm in 1:100
            A_s3[ind_off_s3] = shuffle(A_sur_off)
            p_lv_s3 = (ℵ=A_s3, r=r[sur], N=N_s)
            prob_LV_s3 = ODEProblem(LV_dx!, Ci, tspan, p_lv_s3)
            sol_LV_s3 = solve(prob_LV_s3, AutoVern7(Rodas5()), save_everystep=false, callback=cb)
            bm_LV_s3 = sol_LV_s3.u[length(sol_LV_s3.t)][1:N_s]
            N_LV_s3 = length((1:N_s)[bm_LV_s3.>1.0e-7])

            A_s3_off = [A_s3[i, j] for i in 1:N_s for j in 1:N_s if i != j]
            Jac_s3 = [A_s3[i, j] * C[i] for i in 1:N_s, j in 1:N_s]
            jac_eigen_s3 = eigen(Jac_s3).values
            leading_s3 = jac_eigen_s3[argmax(real.(jac_eigen_s3))]
            append!(all_leading_s3, leading_s3)
            append!(αii_var_s3, var(log.(abs.(diag(A_s3)))))
            append!(αij_var_s3, var(log.(abs.(A_s3_off))))
        end
        next!(progress)
    end
    push!(all_richness_s0_collect_1, rich_s0)
    push!(all_leading_s0_collect_1, all_leading_s0)
    push!(αii_var_s0_collect_1, αii_var_s0)
    push!(αij_var_s0_collect_1, αij_var_s0)
    push!(all_richness_s1_collect_1, rich_s1)
    push!(all_leading_s1_collect_1, all_leading_s1)
    push!(αii_var_s1_collect_1, αii_var_s1)
    push!(αij_var_s1_collect_1, αij_var_s1)
    push!(all_richness_s2_collect_1, rich_s2)
    push!(all_leading_s2_collect_1, all_leading_s2)
    push!(αii_var_s2_collect_1, αii_var_s2)
    push!(αij_var_s2_collect_1, αij_var_s2)
    push!(all_richness_s3_collect_1, rich_s3)
    push!(all_leading_s3_collect_1, all_leading_s3)
    push!(αii_var_s3_collect_1, αii_var_s3)
    push!(αij_var_s3_collect_1, αij_var_s3)
end

@save "../results/randomizing_1.jld2" path_1 all_richness_s0_collect_1 all_leading_s0_collect_1 αii_var_s0_collect_1 αij_var_s0_collect_1 all_richness_s1_collect_1 all_leading_s1_collect_1 αii_var_s1_collect_1 αij_var_s1_collect_1 all_richness_s2_collect_1 all_leading_s2_collect_1 αii_var_s2_collect_1 αij_var_s2_collect_1 all_richness_s3_collect_1 all_leading_s3_collect_1 αii_var_s3_collect_1 αij_var_s3_collect_1
# @save "../results/randomizing_1.jld2" path_1 all_richness_s0_collect_1 all_leading_s0_collect_1 αii_var_s0_collect_1 αij_var_s0_collect_1 all_richness_s1_collect_1 all_leading_s1_collect_1 αii_var_s1_collect_1 αij_var_s1_collect_1 all_richness_s2_collect_1 all_leading_s2_collect_1 αii_var_s2_collect_1 αij_var_s2_collect_1 all_richness_s3_collect_1 all_leading_s3_collect_1 αii_var_s3_collect_1 αij_var_s3_collect_1

#######################################################
########### Analysis from this point on ###############
#######################################################
stability_all_s0 = vcat([sum(real.(all_leading_s0_collect_0[t]) .< 0) / length(path_0) for t in 1:num_temps], [sum(real.(all_leading_s0_collect_1[t]) .< 0) / length(path_1) for t in 1:num_temps])
stability_all_s1 = vcat([sum(real.(all_leading_s1_collect_0[t]) .< 0) / (length(path_0) * 10) for t in 1:num_temps], [sum(real.(all_leading_s1_collect_1[t]) .< 0) / (length(path_1) * 10) for t in 1:num_temps])
stability_all_s2 = vcat([sum(real.(all_leading_s2_collect_0[t]) .< 0) / length(path_0) for t in 1:num_temps], [sum(real.(all_leading_s2_collect_1[t]) .< 0) / length(path_1) for t in 1:num_temps])
stability_all_s3 = vcat([sum(real.(all_leading_s3_collect_0[t]) .< 0) / (length(path_0) * 100) for t in 1:num_temps], [sum(real.(all_leading_s3_collect_1[t]) .< 0) / (length(path_1) * 100) for t in 1:num_temps])

var_ij_all_s0 = vcat([mean(αij_var_s0_collect_0[t]) for t in 1:num_temps], [mean(αij_var_s0_collect_1[t]) for t in 1:num_temps])
var_ij_all_s1 = vcat([mean(αij_var_s1_collect_0[t]) for t in 1:num_temps], [mean(αij_var_s1_collect_1[t]) for t in 1:num_temps])
var_ij_all_s2 = vcat([mean(αij_var_s2_collect_0[t]) for t in 1:num_temps], [mean(αij_var_s2_collect_1[t]) for t in 1:num_temps])
var_ij_all_s3 = vcat([mean(αij_var_s3_collect_0[t]) for t in 1:num_temps], [mean(αij_var_s3_collect_1[t]) for t in 1:num_temps])

var_ii_all_s0 = vcat([mean(αii_var_s0_collect_0[t]) for t in 1:num_temps], [mean(αii_var_s0_collect_1[t]) for t in 1:num_temps])
var_ii_all_s1 = vcat([mean(αii_var_s1_collect_0[t]) for t in 1:num_temps], [mean(αii_var_s1_collect_1[t]) for t in 1:num_temps])
var_ii_all_s2 = vcat([mean(αii_var_s2_collect_0[t]) for t in 1:num_temps], [mean(αii_var_s2_collect_1[t]) for t in 1:num_temps])
var_ii_all_s3 = vcat([mean(αii_var_s3_collect_0[t]) for t in 1:num_temps], [mean(αii_var_s3_collect_1[t]) for t in 1:num_temps])

var_ij_err_s0 = vcat([std(αij_var_s0_collect_0[t]) / sqrt(length(αij_var_s0_collect_0[t])) for t in 1:num_temps], [std(αij_var_s0_collect_1[t]) / sqrt(length(αij_var_s0_collect_1[t])) for t in 1:num_temps])
var_ij_err_s1 = vcat([std(αij_var_s1_collect_0[t]) / sqrt(length(αij_var_s1_collect_0[t])) for t in 1:num_temps], [std(αij_var_s1_collect_1[t]) / sqrt(length(αij_var_s1_collect_1[t])) for t in 1:num_temps])
var_ij_err_s2 = vcat([std(αij_var_s2_collect_0[t]) / sqrt(length(αij_var_s2_collect_0[t])) for t in 1:num_temps], [std(αij_var_s2_collect_1[t]) / sqrt(length(αij_var_s2_collect_1[t])) for t in 1:num_temps])
var_ij_err_s3 = vcat([std(αij_var_s3_collect_0[t]) / sqrt(length(αij_var_s3_collect_0[t])) for t in 1:num_temps], [std(αij_var_s3_collect_1[t]) / sqrt(length(αij_var_s3_collect_1[t])) for t in 1:num_temps])

var_ii_err_s0 = vcat([std(αii_var_s0_collect_0[t]) / sqrt(length(αii_var_s0_collect_0[t])) for t in 1:num_temps], [std(αii_var_s0_collect_1[t]) / sqrt(length(αii_var_s0_collect_1[t])) for t in 1:num_temps])
var_ii_err_s1 = vcat([std(αii_var_s1_collect_0[t]) / sqrt(length(αii_var_s1_collect_0[t])) for t in 1:num_temps], [std(αii_var_s1_collect_1[t]) / sqrt(length(αii_var_s1_collect_1[t])) for t in 1:num_temps])
var_ii_err_s2 = vcat([std(αii_var_s2_collect_0[t]) / sqrt(length(αii_var_s2_collect_0[t])) for t in 1:num_temps], [std(αii_var_s2_collect_1[t]) / sqrt(length(αii_var_s2_collect_1[t])) for t in 1:num_temps])
var_ii_err_s3 = vcat([std(αii_var_s3_collect_0[t]) / sqrt(length(αii_var_s3_collect_0[t])) for t in 1:num_temps], [std(αii_var_s3_collect_1[t]) / sqrt(length(αii_var_s3_collect_1[t])) for t in 1:num_temps])

f = Figure(fontsize=35, size=(1200, 900));
ax = Axis(f[1, 1], xlabel=L"σ(log(α))_{feas}", ylabel="p(Stability)", xlabelsize=50, ylabelsize=50)
scatter!(ax, var_ii, stability_all, color="#FA8328", markersize=15, alpha=0.8, label="αᵢᵢ")
for (x, y, e) in zip(var_ii, stability_all, var_ii_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color=("#FA8328", 0.4), linewidth=1)
    # Horizontal caps
    cap_length_0 = 0.001 * mean(stability_all)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_0, y + cap_length_0], color=("#FA8328", 0.4), linewidth=1)
    lines!(ax, [x + e, x + e], [y - cap_length_0, y + cap_length_0], color=("#FA8328", 0.4), linewidth=1)
end
scatter!(ax, var_ij, stability_all, color="#015845", markersize=15, alpha=0.8, label="αᵢⱼ")
for (x, y, e) in zip(var_ij, stability_all, var_ij_err)
    # Vertical line
    lines!(ax, [x - e, x + e], [y, y], color=("#015845", 0.4), linewidth=1)
    # Horizontal caps
    cap_length_1 = 0.001 * mean(stability_1)  # Length of horizontal caps
    # lines!(ax, [x - cap_length, x + cap_length], [y - e, y - e], color = ("#015845", 0.4), linewidth = 1)
    # lines!(ax, [x - cap_length, x + cap_length], [y + e, y + e], color = ("#015845", 0.4), linewidth = 1)
    lines!(ax, [x - e, x - e], [y - cap_length_1, y + cap_length_1], color=("#015845", 0.4), linewidth=1)
    lines!(ax, [x + e, x + e], [y - cap_length_1, y + cap_length_1], color=("#015845", 0.4), linewidth=1)
end
axislegend(position=:rt)
Label(f[1, 1, TopLeft()], "(b)")
f
