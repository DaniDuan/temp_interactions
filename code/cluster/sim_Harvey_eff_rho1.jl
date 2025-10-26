include("./sim_frame.jl");

N=100
M=50
L = fill(0.3, N)
niche = fill(1.0, M, N)
### Temp params 
num_temps = 31
# ρ_t= [-0.3500 -0.3500]; # realistic covariance
# ρ_t= [0.0000 0.0000]; 
# ρ_t= [-0.9999 -0.9999]; 
ρ_t= [-0.5000 -0.5000]; 
Tr=273.15+10; Ed=3.5 
input_type = ["constant", "leaching", "chemostat", "self-renewing"]

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

# progress = Progress(num_temps; desc="Progress running:")
rich = Float64[]; all_sur = Vector{Vector{Float64}}(); 
all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); all_ℵij_d = Vector{Vector{Float64}}(); all_uℵij = Vector{Vector{Float64}}(); all_lℵij = Vector{Vector{Float64}}();
all_r = Vector{Vector{Float64}}(); 
all_u =  Vector{Vector{Float64}}(); all_m =  Vector{Vector{Float64}}(); 
all_Eu =  Vector{Vector{Float64}}(); all_Em =  Vector{Vector{Float64}}(); 
all_Tpu =  Vector{Vector{Float64}}(); all_Tpm =  Vector{Vector{Float64}}(); 
all_Rrela = Vector{Vector{Float64}}(); all_Crela = Vector{Vector{Float64}}(); all_R =  Vector{Vector{Float64}}(); all_C = Vector{Vector{Float64}}();
all_jac = Vector{Vector{Float64}}(); all_leading = ComplexF64[];

for i in range(0, stop = 30, length = 31)
    try
    T = 273.15 + i
    # next!(progress)

    # N = 1; M= 1
    # Random.seed!(1)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, niche = niche, input_type = input_type[1])

    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    # biomass = [x[1:N] for x in sol.u]
    # resource = [x[N+1:N+M] for x in sol.u]

    sur = (1:N)[bm .> 1.0e-7]
    N_s = length(sur)
    R_t = sol.u[length(sol.t)][N+1:N+M]
    R_rela = R_t ./ sum(R_t)
    C_t = bm[bm .> 1.0e-7]
    C_rela = C_t ./ sum(C_t)
    # CUE
    ϵ = (p.u * x0[N+1:N+M] .* (1 .- p.L) .- p.m) ./ (p.u * x0[N+1:N+M])

    t_len = length(sol.t)
    Rs = sum(x0[N+1:N+M] .+ t_len .- R_t)
    Cs = sum(C_t) .- sum(x0[1:N])
    community_CUE = Cs/Rs
    ## getting effective LV coefficients
    p_lv = Eff_LV_params(p=p, sol=sol);
    r = p_lv.r
    m = p.m
    u = sum(p.u, dims =2)
    # mean E Tp for u and m 
    Eu = p.E[:,1]; Em = p.E[:,2]
    Tpu = p.Tp[:,1]; Tpm = p.Tp[:,2]
    # Resource uptake of survivors 

    ℵii = diag(p_lv.ℵ)
    ℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i != j]
    ℵij_d = [p_lv.ℵ[i, j]/diag(p_lv.ℵ)[i] for i in 1:N for j in 1:N if i != j]
    uℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j > i]
    lℵij = [p_lv.ℵ[j, i] for i in 1:N for j in 1:N if j > i]

    p_lv = Eff_LV_params(p=p, sol=sol);
    LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
    jac_eigen = eigen(LV_jac).values
    jac = vcat(LV_jac...)
    leading = jac_eigen[argmax(real.(jac_eigen))]

    push!(rich, N_s); push!(all_sur, sur); 
    push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_ℵij_d, ℵij_d); push!(all_uℵij, uℵij); push!(all_lℵij, lℵij);
    push!(all_r, r); 
    push!(all_u, vec(u)); push!(all_m, m); 
    push!(all_Eu, Eu); push!(all_Em, Em); 
    push!(all_Tpu, Tpu); push!(all_Tpm, Tpm); 
    push!(all_Rrela, R_rela); push!(all_Crela, C_rela); push!(all_R, R_t); push!(all_C, C_t);
    push!(all_jac, jac); push!(all_leading, leading); 
    catch e 
    push!(rich, NaN); push!(all_sur, NaN); 
    push!(all_ℵii, NaN); push!(all_ℵij, NaN); push!(all_ℵij_d, NaN); push!(all_uℵij, NaN); push!(all_lℵij, NaN);
    push!(all_r, NaN); 
    push!(all_u, NaN); push!(all_m, NaN); 
    push!(all_Eu, NaN); push!(all_Em, NaN); 
    push!(all_Tpu, NaN); push!(all_Tpm, NaN); 
    push!(all_Rrela, NaN); push!(all_Crela, NaN); push!(all_R, NaN); push!(all_C, NaN);
    push!(all_jac, NaN); push!(all_leading, NaN); 
    end 
end 

# R"library(beepr); beep(sound = 4, expr = NULL)"

@save "../data/20250917/p05_test_input/Eff_iters0_ρ1_$(index).jld2" rich all_sur all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_r all_u all_m all_Eu all_Em all_Tpu all_Tpm all_Rrela all_Crela all_R all_C all_jac all_leading
