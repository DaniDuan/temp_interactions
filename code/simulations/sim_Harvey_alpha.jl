include("./sim_frame.jl");

N=100
M=50
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

# Retrieve the environment variable as a string
index_str = ENV["SLURM_ARRAY_TASK_ID"]
# Convert the string to a numeric value (e.g., Integer)
index = parse(Int, index_str)

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
m = fill(1, N)
L = fill(0.3, N)
l = fl(N, M, L)
λ = reshape(sum(l , dims = 3), N, M)
ρ = fill(1, M); ω = fill(0.0, M)
p = (N=N, M=M, u = u, m= m, l = l, ρ =ρ, ω = ω, λ = λ)
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
length(bm .> 1.0e-7)
p_lv = Eff_LV_params(p=p, sol=sol);

all_alpha = vcat(p_lv.ℵ...)
uℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j > i]
lℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j < i]

r = p_lv.r

@save "../results/20240902/uniform03/uniform03_$(index).jld2" all_alpha uℵij lℵij r