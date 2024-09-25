include("./sim_frame.jl");

N=100
M=50
L = fill(0.3, N)
### Temp params 
# num_temps = 31
# Tr=273.15+10; Ed=3.5 
# # ρ_t= [-0.9999 -0.9999]; 
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

# function randtemp_param(N, kw)
#     @unpack T, Tr, var, Ed = kw
#     B0 = [-0.8116 -1.4954]# L = 0.3; mean(CUE0) = 0.22; median(CUE0) = 0.20
#     # B0 = [0.7612 -1.4954]# L = 0.7; mean(CUE0) = 0.22; median(CUE0) = 0.20
#     E_mean = [0.8146 0.5741]; 
#     meanv = [B0 ; E_mean]
#     allu = [rand(Normal(meanv[1,1], var), N),rand(Normal(meanv[2,1], var), N)]
#     allm = [rand(Normal(meanv[1,2], var), N),rand(Normal(meanv[2,2], var), N)]
#     B = [exp.(allu[1]) exp.(allm[1])]
#     E = [allu[2] allm[2]]
    
#     Tpu = 273.15 .+ rand(Normal(35, 5), N)
#     Tpm = Tpu .+ 3
#     Tp = [Tpu Tpm]
#     return B,E,Tp
# end 

ω = fill(0.0, M)
ρ = fill(1, M)
# 0, 0.05, 0.1, 0.2
# T = 273.15 + 10
var_v = range(0, 0.2, 50)
all_sur = Vector{Vector{Float64}}(); all_ℵ = Vector{Vector{Float64}}(); all_r =  Vector{Vector{Float64}}(); all_leading = ComplexF64[]; all_diag_dom = Float64[];
all_C = Vector{Vector{Float64}}(); all_R = Vector{Vector{Float64}}();
for i in 1:length(var_v)
    l = def_l(N, M, L)
    λ = reshape(sum(l , dims = 3), N, M)
    u_sum = rand(Normal(1,  var_v[i]), N)
    u = transpose(rand(Dirichlet(ones(M)),N)).* u_sum
    m = rand(Normal(1,  var_v[i]), N)
    p = (N=N, M=M, l=l, u=u, m=m, ρ = ρ, ω = ω, λ = λ)
    # p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, var = var_v[i])
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol = solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

    p_lv = Eff_LV_params(p=p, sol=sol);

    N = 10
    ℵ = rand(Normal(0, range(0, 1, 50)[i]), N, N)
    ℵ[diagind(ℵ)] .= -1
    r = ones(N)
    p_lv = (ℵ= ℵ, r=r, N=N)
    Ci = fill(0.1, N)
    prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
    sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm_LV = sol_LV.u[length(sol_LV.t)]
    sum(bm_LV .> 1.0e-7)
    
    LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
    jac_eigen = eigen(LV_jac).values
    leading = jac_eigen[argmax(real.(jac_eigen))]
    jac_diag = diag(LV_jac)
    jac_off = [sum(abs.(LV_jac[i, j]) for j in 1:N if j != i) for i in 1:N ]
    diag_dom = sum(abs.(jac_diag) - jac_off .> 0)/N

    ## saving 
    ℵ = vcat(p_lv.ℵ...)
    r = p_lv.r
    C = sol.u[length(sol.t)][1:N]
    sur = (1:N)[C .> 1.0e-7]
    R = sol.u[length(sol.t)][N+1:N+M]

    push!(all_sur, sur); push!(all_ℵ, ℵ); push!(all_r, r); push!(all_leading, leading); push!(all_diag_dom, diag_dom); push!(all_C, C), push!(all_R, R)

end 

@save "../data/20240918/v/v_$(index).jld2" all_sur all_ℵ all_r all_leading all_diag_dom all_C all_R