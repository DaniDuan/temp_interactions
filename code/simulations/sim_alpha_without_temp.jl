include("./sim_frame.jl");
using ProgressMeter, RCall

N=100
M=50
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

# niche_diff = ones(M, N) # 1000 for concentrating, 1.0 for uniform
function fu_niche(N, M, niche_diff)
    u_sum = fill(2.5, M)'
    diri = zeros(Float64, N, M)
    for i in 1:N
        diri[i,:] = rand(Dirichlet(niche_diff[:,i]),1)
    end
    # diri = rand(Dirichlet(ones(M)),N)'
    u = diri.*u_sum
    return u
end

function fl(N, M, L)
    l = zeros(N, M, M)
    ϕ = fill(1.0, M)
    dD = Dirichlet(ϕ[:])
    for i = 1:N
        for α = 1:M
            l[i, α, :] = rand(dD) * L[i]
        end
    end
    return l
end

niche_over = fill(100.0, M, N)
niche_rand = fill(1.0, M, N)
L = fill(0.3, N)

progress = Progress(200; desc="Progress running:")

all_alpha_diff = Float64[]
next!(progress)
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

# u = fu_niche(N, M, niche_diff) # Completely differentiated 
# u = fu_niche(N, M, niche_over) # Complete overlap
u = fu_niche(N, M, niche_rand) # Uniform distribution
m = fill(1, N)
l = fl(N, M, L)
λ = reshape(sum(l , dims = 3), N, M)
ρ = fill(1, M); ω = fill(0.0, M)
p = (N=N, M=M, u = u, m= m, l = l, ρ =ρ, ω = ω, λ = λ)

prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = true, callback=cb)

sol_results = reshape(vcat(sol.u...), N+M, length(sol.u))'
bm = sol_results[:, 1:N]
rbm = sol_results[:, N+1:N+M]
######################### Running the Effective LV
p_lv = Eff_LV_params(p=p, sol=sol);

## running LV
prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
bm_LV = reshape(vcat(sol_LV.u...), N, length(sol_LV.u))'

lines(1:size(bm_LV)[1], bm_LV[:,1])


f = Figure(fontsize = 35, size = (2000, 900));
ax1 = Axis(f[1,1], xlabel = "Time", ylabel = "Consumer Abundance", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,2], xlabel = "Time", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
for i in 1:N 
    lines!(ax1, 1:size(bm)[1], bm[:,i], color = ("#0758AE",0.8), linewidth = 3)
    lines!(ax2, 1:size(bm_LV)[1], bm_LV[:,i], color = ("#4F363E", 0.8), linewidth = 3)
end
l1 = [LineElement(color = ("#0758AE",0.8), linestyle = nothing, linewidth = 3)]
l2 = [LineElement(color = ("#4F363E", 0.8), linestyle = nothing, linewidth = 3)]
Legend(f[1,1], [l1], tellheight = false, tellwidth = false, ["MCM"], halign = :left, valign = :top)
Legend(f[1,2], [l2], tellheight = false, tellwidth = false, ["ELV"], halign = :left, valign = :top)
f
save("../results/MCM_ELV_dynamics.png", f) 





###### 
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
length(bm .> 1.0e-7)
p_lv = Eff_LV_params(p=p, sol=sol);
all_alpha_diff = vcat(p_lv.ℵ...)
hist(all_alpha_diff, bins = 100)
R"library(beepr); beep(sound = 4, expr = NULL)"



######### loading alpha data ############
using Glob
path_over = glob("over*", "../results/over03/")
path_diff = glob("diff*", "../results/diff03/")
path_uniform = glob("uniform*", "../results/uniform03/")

all_alpha_over = Float64[]
progress = Progress(length(path_over); desc="Progress running:")
for index in 1: length(path_over)
    next!(progress)
    @load path_over[index] all_alpha
    append!(all_alpha_over, all_alpha)
end 

all_alpha_diff = Float64[]
progress = Progress(length(path_diff); desc="Progress running:")
for index in 1: length(path_diff)
    next!(progress)
    @load path_diff[index] all_alpha
    append!(all_alpha_diff, all_alpha)
end 

all_alpha_uniform = Float64[]
progress = Progress(length(path_uniform); desc="Progress running:")
for index in 1: length(path_uniform)
    next!(progress)
    @load path_uniform[index] all_alpha
    append!(all_alpha_uniform, all_alpha)
end 

# High niche overlap; High niche differentiation
f = Figure(fontsize = 35, size = (2000, 900));
ax1 = Axis(f[1,1], title = "High niche differentiation", xlabel = "α", ylabel = "frequency", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,2], title = "Uniformly distributed niche", xlabel = "α", ylabel = "frequency", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax3 = Axis(f[1,3], title = "High niche overlap", xlabel = "α", ylabel = "frequency", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
hist!(ax1, all_alpha_diff, color = "#C25E8B", bins = 100)
hist!(ax2, all_alpha_uniform, color = "#82AC6D", bins = 100)
hist!(ax3, all_alpha_over, color = "#5676A5", bins = 100)
f
save("../results/α_niche_03.png", f) 
