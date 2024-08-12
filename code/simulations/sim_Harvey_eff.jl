include("./sim_frame.jl");

N=100
M=50
L = 0.3
### Temp params 
num_temps = 31
ρ_t= [0.0000 0.0000]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

# Retrieve the environment variable as a string
index_str = ENV["SLURM_ARRAY_TASK_ID"]
# Convert the string to a numeric value (e.g., Integer)
index = parse(Int, index_str)

all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); all_ℵij_d = Vector{Vector{Float64}}(); all_uℵij = Vector{Vector{Float64}}(); all_lℵij = Vector{Vector{Float64}}();
all_ℵii_sur =  Vector{Vector{Float64}}(); all_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); all_ℵij_d_sur =  Vector{Vector{Union{Float64, Missing}}}(); all_uℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); all_lℵij_sur = Vector{Vector{Union{Float64, Missing}}}();
all_r = Vector{Vector{Float64}}(); 
all_leading = Float64[]; all_diag = Vector{Vector{Float64}}();radi = Vector{Vector{Float64}}(); diag_dominance = Float64[];
all_u =  Vector{Vector{Float64}}(); all_m =  Vector{Vector{Float64}}(); RO = Vector{Vector{Union{Float64, Missing}}}(); ulO = Vector{Vector{Union{Float64, Missing}}}(); Rul = Vector{Vector{Union{Float64, Missing}}}(); 
all_Eu =  Vector{Vector{Float64}}(); all_Em =  Vector{Vector{Float64}}(); all_Eu_sur = Vector{Vector{Float64}}(); all_Em_sur = Vector{Vector{Float64}}();
all_Tpu =  Vector{Vector{Float64}}(); all_Tpm =  Vector{Vector{Float64}}(); all_Tpu_sur = Vector{Vector{Float64}}(); all_Tpm_sur = Vector{Vector{Float64}}()

@time for i in range(0, stop = 30, length = 31)
    T = 273.15 + i

    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    sur = (1:N)[bm .> 1.0e-7]
    N_s = length(sur)
    R_t = sol.u[length(sol.t)][N+1:N+M]
    C_t = sol.u[length(sol.t)][1:N][sur]
    ## getting effective LV coefficients
    p_lv = Eff_LV_params(p=p, sol=sol);
    # number of species with r>0 at equilibium 
    N_sur = sum(p_lv.r .> 0)
    sur_r = p_lv.r[p_lv.r.>0]
    # mean uptake and respiration 
    m = p.m[p_lv.r.>0]
    u = sum(p.u, dims =2)[p_lv.r.>0]
    # mean E Tp for u and m 
    Eu = p.E[:,1]; Em = p.E[:,2]
    Eu_sur = Eu[sur]; Em_sur = Em[sur]
    Tpu = p.Tp[:,1]; Tpm = p.Tp[:,2]
    Tpu_sur = Tpu[sur]; Tpm_sur = Tpm[sur]
    # Resource uptake of survivors 
    u_sur = p.u[sur,:]

    # survivor interaction coefficients
    # ℵ = p_lv.ℵ # interactions
    sur_ℵ = p_lv.ℵ[sur, sur] # interactions in the possibily surviving community
    ℵii = diag(p_lv.ℵ)
    ℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if i != j]
    ℵij_d = [p_lv.ℵ[i, j]/diag(p_lv.ℵ)[i] for i in 1:N for j in 1:N if i != j]
    uℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j > i]
    lℵij = [p_lv.ℵ[i, j] for i in 1:N for j in 1:N if j < i]
    # eigenvalue for jacobian 
    LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
    jac_eigen = eigen(LV_jac).values
    leading = jac_eigen[argmax(real.(jac_eigen))]
    jac_diag = diag(LV_jac)
    jac_off = [sum(abs.(LV_jac[i, j]) for j in 1:N if j != i) for i in 1:N ]
    diag_dom = sum(abs.(jac_diag) - jac_off .> 0)/N
    if N_s > 1
        u_tR = mapslices(x -> x .* R_t, u_sur, dims=2) # getting the actual uptake
        u_t = mapslices(x -> x .* C_t, u_tR, dims=1) # getting the actual uptake
        R_over = 1 .- [bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
        l_t = p.l[sur,:,:]
        ul = zeros(Float64, N_s, M)
        for s in 1: N_s
                uli = zeros(Float64, M, M)
                for α in 1:M
                uli[α,:] = u_t[s, α] .* l_t[s, α, :]
                end 
                ul[s,:] = sum(uli, dims = 1)
        end 
        ul_over = 1 .- [bray_curtis_dissimilarity(ul[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
        Rul_over = ul_over - R_over
        ℵij_sur = [sur_ℵ[i, j] for i in 1:N_s for j in 1:N_s if i != j]
        ℵij_d_sur = [sur_ℵ[i, j]/diag(sur_ℵ)[i] for i in 1:N_s for j in 1:N_s if i != j]
        uℵij_sur = [sur_ℵ[i, j] for i in 1:N_s for j in 1:N_s if j > i]
        lℵij_sur = [sur_ℵ[i, j] for i in 1:N_s for j in 1:N_s if j < i]
        # UDLD_sur = [sur_ℵ[i, j]/sur_ℵ[j, i] for i in 1:N_s for j in 1:N_s if j != i]
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_ℵij_d, ℵij_d); push!(all_uℵij, uℵij); push!(all_lℵij, lℵij);
        push!(all_ℵii_sur, diag(sur_ℵ)); push!(all_ℵij_sur, ℵij_sur); push!(all_ℵij_d_sur, ℵij_d_sur); push!(all_uℵij_sur, uℵij_sur); push!(all_lℵij_sur, lℵij_sur); 
        push!(all_r, sur_r);
        push!(all_u, u); push!(all_m, m); push!(RO, R_over); push!(ulO, ul_over); push!(Rul, Rul_over);
        push!(all_Eu, Eu); push!(all_Em, Em); push!(all_Eu_sur, Eu_sur); push!(all_Em_sur, Em_sur);
        push!(all_Tpu, Tpu); push!(all_Tpm, Tpm); push!(all_Tpu_sur, Tpu_sur); push!(all_Tpm_sur, Tpm_sur);
        push!(all_leading, leading); push!(all_diag, jac_diag); push!(radi, jac_off); push!(diag_dominance, diag_dom)
    else 
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_ℵij_d, ℵij_d); push!(all_uℵij, uℵij); push!(all_lℵij, lℵij);
        push!(all_ℵii_sur, diag(sur_ℵ)); push!(all_ℵij_sur, [missing]); push!(all_ℵij_d_sur, [missing]); push!(all_uℵij_sur, [missing]); push!(all_lℵij_sur, [missing]); 
        push!(all_r, sur_r);
        push!(all_u, u); push!(all_m, m); push!(RO, [missing]); push!(ulO, [missing]); push!(Rul, [missing]);
        push!(all_Eu, Eu); push!(all_Em, Em); push!(all_Eu_sur, Eu_sur); push!(all_Em_sur, Em_sur);
        push!(all_Tpu, Tpu); push!(all_Tpm, Tpm); push!(all_Tpu_sur, Tpu_sur); push!(all_Tpm_sur, Tpm_sur);
        push!(all_leading, leading); push!(all_diag, jac_diag); push!(radi, jac_off); push!(diag_dominance, diag_dom)
    end
end 

@save "../data/20240809/p0/Eff_iters0_$(index).jld2" all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_ℵii_sur all_ℵij_sur all_ℵij_d_sur all_uℵij_sur all_lℵij_sur all_r  all_u all_m RO ulO Rul all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_diag radi diag_dominance 

# @load "../data/Eff_iter/Eff_iters1.jld2" all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_ℵii_sur all_ℵij_sur all_ℵij_d_sur all_uℵij_sur all_lℵij_sur all_r  all_u all_m RO ulO Rul all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_diag radi diag_dominance 
