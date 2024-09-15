include("./sim_frame.jl");
using ProgressMeter, RCall

N=50
M=50
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

function fu_niche(N, M, niche)
    u_sum = fill(2.5, N)
    # u_sum = rand(truncated(Normal(1.58, 1), 0, Inf), N)'
    # maximum(u_sum)
    diri = zeros(Float64, N, M)
    for i in 1:N
        diri[i,:] = rand(Dirichlet(niche[:,i]),1)
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
# u = fill(2.5/M, N, M)
m = fill(1, N)
L = fill(0.3, N)
l = fl(N, M, L)
λ = reshape(sum(l , dims = 3), N, M)
ρ = fill(1, M); ω = fill(0.0, M)
p = (N=N, M=M, u = u, m= m, l = l, ρ =ρ, ω = ω, λ = λ)
prob = ODEProblem(dxx!, x0, tspan, p)
# sol =solve(prob, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
# bm = sol.u[length(sol.t)][:, 1:N]
# rbm = sol.u[length(sol.t)][:, N+1:N+M]

sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
rbm = sol.u[length(sol.t)][N+1:N+M]

######################### Running the Effective LV
p_lv = Eff_LV_params(p=p, sol=sol);

## running LV
prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
# sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
# bm_LV = reshape(vcat(sol_LV.u...), N, length(sol_LV.u))'

sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm_LV = sol_LV.u[length(sol_LV.t)]

all_alpha = vcat(p_lv.ℵ...)
uℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j > i]
lℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j < i]

r = p_lv.r

# lines(1:size(bm_LV)[1], bm_LV[:,1])


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
save("../results/MCM_ELV_dynamics.pdf", f) 


################## heatmap ######################
N = 100; M = 50
niche_over = fill(100.0, M, N)
niche_rand = fill(1.0, M, N)

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
u_over = fu_niche(N, M, niche_over) # Complete overlap
u_rand = fu_niche(N, M, niche_rand) # Uniform distribution


# combined_data = [u_rand; u_over]  
# color_limits = extrema(combined_data)
f = Figure(fontsize = 35, size = (2100, 700));
ax1 = Axis(f[1,1][1,1], title = "High niche differentiation", ygridvisible = false, xgridvisible = false, xlabelsize = 35, ylabelsize = 35)
ax2 = Axis(f[1,2][1,1], title = "Uniformly distributed niche", ygridvisible = false, xgridvisible = false, xlabelsize = 35, ylabelsize = 35)
ax3 = Axis(f[1,3][1,1], title = "High niche overlap", ygridvisible = false, xgridvisible = false, xlabelsize = 35, ylabelsize = 35)
hm1 = heatmap!(ax1, u_diff, colormap = Reverse(:grayC))
hm2 = heatmap!(ax2, u_rand, colormap = Reverse(:grayC)) # , colorrange = color_limits
hm3 = heatmap!(ax3, u_over, colormap = Reverse(:grayC))
ax1.xticks = ([], []); ax1.yticks = ([], [])
ax2.xticks = ([], []); ax2.yticks = ([], [])
ax3.xticks = ([], []); ax3.yticks = ([], [])
Colorbar(f[1,1][1,2], hm1)  
Colorbar(f[1,2][1,2], hm2)  
Colorbar(f[1,3][1,2], hm3)  
f
save("../results/α_niche_u03.pdf", f) 


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

all_alpha_over = Float64[]; all_αii_over = Float64[]; all_αij_over = Float64[];
all_upper_over = Float64[]; all_lower_over = Float64[];
progress = Progress(length(path_over); desc="Progress running:")
for index in 1: length(path_over)
    next!(progress)
    @load path_over[index] all_alpha
    A = reshape(all_alpha, N, N)
    αii = diag(A)
    αij = [A[i, j] for i in 1:N for j in 1:N if i != j]
    append!(all_alpha_over, all_alpha); append!(all_αii_over, αii);  append!(all_αij_over, αij);
    append!(all_upper_over, uℵij); append!(all_lower_over, lℵij)
end 

all_alpha_diff = Float64[]; all_αii_diff = Float64[]; all_αij_diff = Float64[];
all_upper_diff = Float64[]; all_lower_diff = Float64[];
all_r_diff = Float64[]
progress = Progress(length(path_diff); desc="Progress running:")
for index in 1: length(path_diff)
    next!(progress)
    @load path_diff[index] all_alpha
    A = reshape(all_alpha, N, N)
    αii = diag(A)
    αij = [A[i, j] for i in 1:N for j in 1:N if i != j]
    append!(all_alpha_diff, all_alpha); append!(all_αii_diff, αii);  append!(all_αij_diff, αij);
    append!(all_upper_diff, uℵij); append!(all_lower_diff, lℵij); 

end 

all_alpha_uniform = Float64[]; all_αii_uniform = Float64[]; all_αij_uniform = Float64[];
all_upper_uniform = Float64[]; all_lower_uniform = Float64[];
progress = Progress(length(path_uniform); desc="Progress running:")
for index in 1: length(path_uniform)
    next!(progress)
    @load path_uniform[index] all_alpha uℵij lℵij
    A = reshape(all_alpha, N, N)
    αii = diag(A)
    αij = [A[i, j] for i in 1:N for j in 1:N if i != j]
    append!(all_alpha_uniform, all_alpha); append!(all_αii_uniform, αii);  append!(all_αij_uniform, αij);
    append!(all_upper_uniform, uℵij); append!(all_lower_uniform, lℵij)
end 

# High niche overlap; High niche differentiation #82AC6D; #C25E8B
f = Figure(fontsize = 35, size = (2100, 900));
ax1 = Axis(f[1,1], title = "High niche differentiation", xlabel = "α", ylabel = "frequency (αii)", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax2)
hist!(ax1, all_αii_diff, color = ("#FA8328", 0.7), bins = 100)
hist!(ax2, all_αij_diff, color = ("#069F66", 0.7), bins = 100)
linkxaxes!(ax1,ax2)
ax3 = Axis(f[1,2], title = "Uniformly distributed niche", xlabel = "α", ylabel = "", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax4 = Axis(f[1,2], ylabel = "", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax4)
hist!(ax3, all_αii_uniform, color = ("#FA8328", 0.7), bins = 100)
hist!(ax4, all_αij_uniform, color = ("#069F66", 0.7), bins = 100)
linkxaxes!(ax3,ax4)
ax5 = Axis(f[1,3], title = "High niche overlap", xlabel = "α", ylabel = "", xlabelsize = 35, ylabelsize = 35, ygridvisible = true, xgridvisible = true)
ax6 = Axis(f[1,3], ylabel = "frequency (αij)", xlabelsize = 35, ylabelsize = 35, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hidespines!(ax6)
hist!(ax5, all_αii_over, color = ("#FA8328", 0.5), bins = 100)
hist!(ax6, all_αij_over, color = ("#069F66", 0.5), bins = 100)
linkxaxes!(ax5,ax6)
p1 = [PolyElement(color = ("#FA8328", 0.7), strokecolor = :transparent)]
p2 = [PolyElement(color = ("#069F66", 0.7), strokecolor = :transparent)]
Legend(f[1,1], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
Legend(f[1,2], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
Legend(f[1,3], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
f
save("../results/α_niche_03.pdf", f) 


##########################
f = Figure(fontsize = 35, size = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "α", ylabel = "frequency (αii)", xlabelsize = 50, ylabelsize = 50, ygridvisible = true, xgridvisible = true)
ax2 = Axis(f[1,1], ylabel = "frequency (αij)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), ygridvisible = false, xgridvisible = false, xticklabelsvisible = false, xlabelvisible = false)
hist!(ax1, all_αii_uniform, color = ("#FA8328", 0.7), bins = 100)
hist!(ax2, all_αij_uniform, color = ("#069F66", 0.7), bins = 100)
linkxaxes!(ax1,ax2)
Legend(f[1,1], [p1, p2], tellheight = false, tellwidth = false, ["αii", "αij"], halign = :left, valign = :top)
f
save("../results/α_uniform_03.png", f) 
