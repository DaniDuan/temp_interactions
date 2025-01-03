include("./sim_frame.jl");

N=100
M=50
L = fill(0.3, N)
var_B = range(0, 0.3, 50)
niche = fill(1.0, M, N)

###################################
# Generate MiCRM parameters
tspan = (0.0, 2.5e10)
x0 = vcat(fill(0.1, N), fill(1, M)) 
# here we define a callback that terminates integration as soon as system reaches steady state
# condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
condition(du, t, integrator) = norm(integrator(t, Val{1})[N:N+M]) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

# Retrieve the environment variable as a string
index_str = ENV["SLURM_ARRAY_TASK_ID"]
# Convert the string to a numeric value (e.g., Integer)
index = parse(Int, index_str)

function generate_var_param(N, M, L, var_Bv, niche)
    diri = zeros(Float64, N, M)
    for i in 1:N
        diri[i,:] = rand(Dirichlet(niche[:,i]),1)
    end

    u_sum = rand(Truncated(Normal(2.5, var_Bv), 0, Inf))
    u = diri.*u_sum
    m = rand(Normal(1, var_Bv), N)
    l = def_l(N, M, L)
    ρ = ones(M); ω = zeros(M)
    λ = repeat(L, 1, M)

    return (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ)
end 
    
# 0, 0.05, 0.1, 0.2
all_sur = Vector{Vector{Float64}}(); all_ℵ = Vector{Vector{Float64}}(); all_jac = Vector{Vector{Float64}}(); all_r =  Vector{Vector{Float64}}(); all_leading = ComplexF64[];# all_diag_dom = Float64[];
all_C = Vector{Vector{Float64}}(); all_R = Vector{Vector{Float64}}();
for i in 1:length(var_B)
    # p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, var_B = var_B[i])
    p = generate_var_param(N, M, L, var_B[i], niche)

    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

    p_lv = Eff_LV_params(p=p, sol=sol);
    LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
    jac_eigen = eigen(LV_jac).values
    leading = jac_eigen[argmax(real.(jac_eigen))]
    # jac_diag = diag(LV_jac)
    # jac_off = [sum(abs.(LV_jac[i, j]) for j in 1:N if j != i) for i in 1:N ]
    # diag_dom = sum(abs.(jac_diag) - jac_off .> 0)/N

    ## saving 
    ℵ = vcat(p_lv.ℵ...)
    jac = vcat(LV_jac...)
    r = p_lv.r
    C = sol.u[length(sol.t)][1:N]
    sur = (1:N)[C .> 1.0e-7]
    R = sol.u[length(sol.t)][N+1:N+M]

    push!(all_sur, sur); push!(all_ℵ, ℵ); push!(all_r, r); push!(all_jac, jac); push!(all_leading, leading); push!(all_C, C), push!(all_R, R)

end 

@save "../data/20241218/v/v_$(index).jld2" all_sur all_ℵ all_r all_jac all_leading all_C all_R