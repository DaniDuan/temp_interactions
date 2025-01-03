include("./sim_frame.jl");

N=100
M=50
L = fill(0.3, N)
### Temp params 
num_temps = 31
ρ_t= [-0.3500 -0.3500]; # realistic covariance
# ρ_t= [0.0000 0.0000]; 
# ρ_t= [-0.9999 -0.9999]; 
Tr=273.15+10; Ed=3.5 
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

### Testing L
# l_v = 0.1 * ((index-1)%7+1)
# L = fill(l_v, N)
### Testing niche
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
niche_all =(niche_over, niche_rand, niche_diff)
niche =niche_rand #niche_all[(index-1)%3+1]

# progress = Progress(num_temps; desc="Progress running:")
rich = Float64[]; all_sur = Vector{Vector{Float64}}(); all_ϵ =  Vector{Vector{Float64}}();
all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); all_ℵij_d = Vector{Vector{Float64}}(); all_uℵij = Vector{Vector{Float64}}(); all_lℵij = Vector{Vector{Float64}}();
all_ℵii_sur =  Vector{Vector{Float64}}(); all_ℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); all_ℵij_d_sur =  Vector{Vector{Union{Float64, Missing}}}(); all_uℵij_sur = Vector{Vector{Union{Float64, Missing}}}(); all_lℵij_sur = Vector{Vector{Union{Float64, Missing}}}();
all_r = Vector{Vector{Float64}}(); all_r_sur = Vector{Vector{Float64}}();
all_leading = ComplexF64[]; all_H_leading = ComplexF64[]; all_Jac = Vector{Vector{ComplexF64}}(); diag_dominance = Float64[];
all_u =  Vector{Vector{Float64}}(); all_m =  Vector{Vector{Float64}}(); all_u_sur = Vector{Vector{Float64}}()
# RO =  Vector{Vector{Float64}}(); ulO =  Vector{Vector{Float64}}(); Rul =  Vector{Vector{Float64}}(); 
# RO_sur = Vector{Vector{Union{Float64, Missing}}}(); ulO_sur = Vector{Vector{Union{Float64, Missing}}}(); Rul_sur = Vector{Vector{Union{Float64, Missing}}}(); 
all_Eu =  Vector{Vector{Float64}}(); all_Em =  Vector{Vector{Float64}}(); all_Eu_sur = Vector{Vector{Float64}}(); all_Em_sur = Vector{Vector{Float64}}();
all_Tpu =  Vector{Vector{Float64}}(); all_Tpm =  Vector{Vector{Float64}}(); all_Tpu_sur = Vector{Vector{Float64}}(); all_Tpm_sur = Vector{Vector{Float64}}();
all_Rrela = Vector{Vector{Float64}}(); all_Crela = Vector{Vector{Float64}}(); all_R =  Vector{Vector{Float64}}(); all_C = Vector{Vector{Float64}}();
all_com_CUE = Float64[]

for i in range(0, stop = 30, length = 31)
    T = 273.15 + i
    # next!(progress)

    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed, niche = niche)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = true, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
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
    r_sur = r[sur]
    # # number of species with r>0 at equilibium 
    # N_sur = sum(p_lv.r .> 0)
    # mean uptake and respiration 
    m = p.m
    u = sum(p.u, dims =2)
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
    lℵij = [p_lv.ℵ[j, i] for i in 1:N for j in 1:N if j > i]
    # # Resource uptake and cross-feeding 
    # u_all = mapslices(x -> x .* x0[1:N], p.u, dims=1) # getting the actual uptake
    # u_all_v = p.u * p.u'
    # R_over = [cosine_similarity(u_all[i,:], u_all[j,:]) for i in 1:N for j in 1:N if j != i] #; mean(R_over)
    # # 1 .- [bray_curtis_dissimilarity(u_all[i,:], u_all[j,:]) for i in 1:N for j in 1:N if j != i]
    # ul = zeros(Float64, N, M)
    # for s in 1: N
    #     uli = zeros(Float64, M, M)
    #     for α in 1:M
    #     uli[α,:] = u_all[s, α] .* p.l[s, α, :]
    #     end 
    #     ul[s,:] = sum(uli, dims = 1)
    # end 
    # ul_over = [cosine_similarity(ul[i,:], u_all[j,:]) for i in 1:N for j in 1:N if j != i] # ; mean(ul_over)
    # ul_over_v = ul * p.u'
    # # 1 .- [bray_curtis_dissimilarity(u_all[i,:], ul[j,:]) for i in 1:N for j in 1:N if j != i]
    # Rul_over = ul_over - u_all

    #### eigenvalue for jacobian 
    LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
    jac_eigen = eigen(LV_jac).values
    leading = jac_eigen[argmax(real.(jac_eigen))]

    ### reactivity # https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2014.00021/full
    LV_H = (LV_jac + LV_jac')./2
    H_eigen = eigen(LV_H).values
    H_leading = H_eigen[argmax(real.(H_eigen))]

    if N_s > 1
        # u_tR = mapslices(x -> x .* R_t, u_sur, dims=2) # getting the actual uptake
        # u_t = mapslices(x -> x .* C_t, u_tR, dims=1) # getting the actual uptake
        # R_over_sur = 1 .- [bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
        # l_t = p.l[sur,:,:]
        # ul_sur = zeros(Float64, N_s, M)
        # for s in 1: N_s
        #         uli = zeros(Float64, M, M)
        #         for α in 1:M
        #         uli[α,:] = u_t[s, α] .* l_t[s, α, :]
        #         end 
        #         ul_sur[s,:] = sum(uli, dims = 1)
        # end 
        # ul_over_sur = 1 .- [bray_curtis_dissimilarity(u_t[i,:], ul_sur[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
        # Rul_over_sur = ul_over_sur - R_over_sur
        ℵij_sur = [sur_ℵ[i, j] for i in 1:N_s for j in 1:N_s if i != j]
        ℵij_d_sur = [sur_ℵ[i, j]/diag(sur_ℵ)[i] for i in 1:N_s for j in 1:N_s if i != j]
        uℵij_sur = [sur_ℵ[i, j] for i in 1:N_s for j in 1:N_s if j > i]
        lℵij_sur = [sur_ℵ[i, j] for i in 1:N_s for j in 1:N_s if j < i]
        # # Jacobian_diag_dom
        jac_diag = diag(LV_jac)
        jac_off = [sum(abs.(LV_jac[i, j]) for j in 1:N_s if j != i) for i in 1:N_s ]
        diag_dom = sum(abs.(jac_diag) - jac_off .> 0)/N_s

        # UDLD_sur = [sur_ℵ[i, j]/sur_ℵ[j, i] for i in 1:N_s for j in 1:N_s if j != i]
        push!(rich, N_s); push!(all_sur, sur); push!(all_ϵ, ϵ);
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_ℵij_d, ℵij_d); push!(all_uℵij, uℵij); push!(all_lℵij, lℵij);
        push!(all_ℵii_sur, diag(sur_ℵ)); push!(all_ℵij_sur, ℵij_sur); push!(all_ℵij_d_sur, ℵij_d_sur); push!(all_uℵij_sur, uℵij_sur); push!(all_lℵij_sur, lℵij_sur); 
        push!(all_r, r); push!(all_r_sur, r_sur);
        push!(all_u, vec(u)); push!(all_m, m); push!(all_u_sur, vec(p.u[sur,1:M])) # reshape(u_sur, 3, 50)
        # push!(RO, R_over); push!(ulO, ul_over); push!(Rul, Rul_over);
        # push!(RO_sur, R_over_sur); push!(ulO_sur, ul_over_sur); push!(Rul_sur, Rul_over_sur);
        push!(all_Eu, Eu); push!(all_Em, Em); push!(all_Eu_sur, Eu_sur); push!(all_Em_sur, Em_sur);
        push!(all_Tpu, Tpu); push!(all_Tpm, Tpm); push!(all_Tpu_sur, Tpu_sur); push!(all_Tpm_sur, Tpm_sur);
        push!(all_leading, leading); push!(all_H_leading, H_leading); push!(all_Jac, vec(LV_jac)); push!(diag_dominance, diag_dom);
        push!(all_Rrela, R_rela); push!(all_Crela, C_rela); push!(all_R, R_t); push!(all_C, C_t);
        push!(all_com_CUE, community_CUE)
    else 
        push!(rich, N_s); push!(all_sur, sur); push!(all_ϵ, ϵ);
        push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_ℵij_d, ℵij_d); push!(all_uℵij, uℵij); push!(all_lℵij, lℵij);
        push!(all_ℵii_sur, diag(sur_ℵ)); push!(all_ℵij_sur, [missing]); push!(all_ℵij_d_sur, [missing]); push!(all_uℵij_sur, [missing]); push!(all_lℵij_sur, [missing]); 
        push!(all_r, all_r); push!(all_r_sur, r_sur);
        push!(all_u, vec(u)); push!(all_m, m); push!(all_u_sur, vec(p.u[sur,1:M])) # reshape(u_sur, 3, 50)
        # push!(RO, [missing]); push!(ulO, [missing]); push!(Rul, [missing]);
        # push!(RO_sur, [missing]); push!(ulO_sur, [missing]); push!(Rul_sur, [missing]);
        push!(all_Eu, Eu); push!(all_Em, Em); push!(all_Eu_sur, Eu_sur); push!(all_Em_sur, Em_sur);
        push!(all_Tpu, Tpu); push!(all_Tpm, Tpm); push!(all_Tpu_sur, Tpu_sur); push!(all_Tpm_sur, Tpm_sur);
        push!(all_leading, leading); push!(all_H_leading, H_leading); push!(all_Jac, vec(LV_jac)); push!(diag_dominance, [missing]);
        push!(all_Rrela, R_rela); push!(all_Crela, C_rela); push!(all_R, R_t); push!(all_C, C_t);
        push!(all_com_CUE, community_CUE)
    end
end 

# R"library(beepr); beep(sound = 4, expr = NULL)"

@save "../data/20250103/p_re/Eff_iters_re_$(index)_n$((index-1)%3+1).jld2" rich all_sur all_ϵ all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_ℵii_sur all_ℵij_sur all_ℵij_d_sur all_uℵij_sur all_lℵij_sur all_r all_r_sur all_u all_m all_u_sur all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_H_leading all_Jac diag_dominance all_Rrela all_Crela all_R all_C all_com_CUE

# @load "../data/20240915/p_re/Eff_iters_re_$(index).jld2" all_ℵii all_ℵij all_ℵij_d all_uℵij all_lℵij all_ℵii_sur all_ℵij_sur all_ℵij_d_sur all_uℵij_sur all_lℵij_sur all_r all_r_sur all_u all_m RO ulO Rul RO_sur ulO_sur Rul_sur all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur all_leading all_H_leading all_diag radi diag_dominance all_Rrela all_Crela all_R all_C all_com_CUE
