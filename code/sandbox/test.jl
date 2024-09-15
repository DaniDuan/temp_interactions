
N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]

tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

T = 273.15 + 25
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
p_lv = Eff_LV_params(p=p, sol=sol);


Ci = fill(0.1, N)
prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
C = sol_LV.u[length(sol_LV.t)][1:N]

LV_Jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
bm = sol.u[length(sol.t)][1:p_lv.N]
sur = (1:p_lv.N)[bm .> 1.0e-7]

Jac_off = [sum(abs.(LV_Jac[i, j]) for j in 1:N if j != i) for i in 1:N ]
diag(LV_Jac) - Jac_off

sum(diag(LV_Jac))
sum(eigen(LV_Jac).values)
sum(diag(LV_Jac) - Jac_off .> 0)

sum(eigen(LV_Jac).values .> 0)
Jsym = 1/2 *(LV_Jac + LV_Jac')
sum(eigen(Jsym).values .> 0)
Jskew = 1/2 *(LV_Jac - LV_Jac')
sum(real(eigen(Jskew).values) .< 1.0e-7)
println(imag(eigen(Jskew).values))

eigen(Jsym).vectors[100,:]
bm = sol.u[length(sol.t)][1:p_lv.N]
sur = (1:p_lv.N)[bm .> 1.0e-7]
N_s = length(sur)

u_sur = p.u[sur,:]
R_t = sol.u[length(sol.t)][N+1:N+M]
C_t = sol.u[length(sol.t)][1:N][sur]
u_tR = mapslices(x -> x .* R_t, u_sur, dims=2) # getting the actual uptake
u_t = mapslices(x -> x .* C_t, u_tR, dims=1) # getting the actual uptake

R_over = 1 .-[bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s if j != i]
R_over_matrix = reshape(1 .-[bray_curtis_dissimilarity(u_t[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s], N_s, N_s)

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
ul_over_matrix = transpose(reshape(1 .-[bray_curtis_dissimilarity(ul[i,:], u_t[j,:]) for i in 1:N_s for j in 1:N_s], N_s, N_s))

est = ul_over - R_over
Plots.histogram(ul_over./R_over)
Plots.histogram(R_over)
Plots.histogram(ul_over)
Plots.histogram(ul_over./R_over)
# sur_ℵ = reshape(sur_ℵ, N_sur*N_sur, 1)
# Plots.scatter(est_ℵ, sur_ℵ)
Plots.histogram(est)


#alpha diversity - Shannon
N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5 #[-0.1384 -0.1384]

tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

T = 273.15 + 25
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol.u[length(sol.t)][1:N]
sur = (1:p_lv.N)[bm .> 1.0e-7]

p = bm[sur]./sum(bm[sur])
Shannon = - sum(p .* log.(p))
Simpson = 1/ sum(p .^2)


#### Script for secondary derivative
using Dierckx

x = x_t
y = log.(abs.(Eff_results.sum_αij))
spl = Spline1D(x, y, k=3)  # k=3 spline interpolation

d2 = [derivative(spl, xi, 2) for xi in x]
plot(d2)

data_ii = DataFrame(y = log.(abs.(Eff_results.αii)), x = x);
B, E= coef(lm(@formula(y ~ x), data_ii))
B = exp(B)
yii = log.(B * exp.(E .* x_t))

f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "α", xlabelsize = 50, ylabelsize = 50)
lines!(ax, x, yii, color = ("#285C93",1), linewidth = 5, label = "")
lines!(ax, x, y, color = ("#E17542",1), linewidth = 5, label = "")
# axislegend(position = :rb)
f

######### Script for finding arrhenius stopping point
# Temp_rich
all_yij = [mean(αij[i]) for i in 1:num_temps]
dataii = DataFrame(y = log.(abs.(all_yij)), x = x_t);
fit1 = lm(@formula(y ~ x), dataii)

X = hcat(ones(num_temps), x_t, x_t.^2)
fit2 = lm(X, log.(abs.(all_yij)))

for i in 1: num_temps
    all_yij = [mean(αij[i]) for i in 1:num_temps]
    t_points = 1: num_temps-i
    dataii = DataFrame(y = log.(abs.(all_yij))[t_points], x = x_t[t_points]);
    fit1 = lm(@formula(y ~ x), dataii)
    X = hcat(ones(length(t_points)), x_t[t_points], x_t[t_points].^2)
    fit2 = lm(X, log.(abs.(all_yij[t_points])))
    if AIC(fit2, length(t_points)) <= AIC(fit1, length(t_points))
        print(i)
        break 
    end 
end 

plot(Temp_rich, log.(abs.(all_yij)))


#################################
c = Int64[]
for j in 1: 9
    for n in 1:10
        i = 10*(j-1) + n
        push!(c, i)
    end 
end 
println(c)

############ Plotting different alpha for Dirichlet ########
# Function to sample from the Dirichlet distribution
function sample_dirichlet(α, num_samples)
    d = Dirichlet(α)
    rand(d, num_samples)
end

alphas = [[100.0, 100.0, 100.0], [1.0, 1.0, 100.0], [1.0, 1.0, 1.0], [100.0, 100.0, 1.0]]
num_samples = 1000

f1 = Figure(size = (800, 800));
f2 = Figure(size = (800, 800));

for i in 1:length(alphas)
    ax1 = Axis(f1[div(i-1, 2)+1, (i-1)%2+1], title = "α = " * string(alphas[i]));
    ax2 = Axis(f2[div(i-1, 2)+1, (i-1)%2+1], title = "α = " * string(alphas[i]));
    samples = sample_dirichlet(alphas[i], num_samples)
    xlims!(ax1, -0.1, 1.1);
    ylims!(ax1, -0.1, 1.1);
    xlims!(ax2, -0.1, 1.1);
    scatter!(ax1, samples[1, :], samples[2, :], color = samples[3, :], colormap = :viridis);
    hist!(ax2, samples[1, :])
    hist!(ax2, samples[2, :])
    hist!(ax2, samples[3, :])
end

display(f1)
display(f2)


############# Define the color gradient from #376298 to #9A2B1A
using Colors

color_start = RGBA(parse(Colorant, "#376298"), 0.4)
color_end = RGBA(parse(Colorant, "#9A2B1A"), 0.4)
color_gradient = range(color_start, stop=color_end, length=length(all_ii_collect))


A = [1 2 3 4; 4 5 6 6; 7 8 9 9; 1 3 3 5]

idx = CartesianIndices(A)
ind_list = collect(idx)
ind_off = [ind_list[i,j] for i in 1:N for j in 1:N if i != j]
A[ind_off[1]]



##############
progress = Progress(20* 6; desc="Progress running:")

leading_eigens =  Vector{Vector{ComplexF64}}()
for j in 1:20
    all_leading = ComplexF64[]
    for i in range(0, stop = 30, length = 6)
        T = 273.15 + i
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        ## run simulation
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        p_lv = Eff_LV_params(p=p, sol=sol);
        LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
        jac_eigen = eigen(LV_jac).values
        leading = jac_eigen[argmax(real.(jac_eigen))]
        push!(all_leading, leading)
        next!(progress)
    end
    push!(leading_eigens, all_leading)
end
R"library(beepr); beep(sound = 4, expr = NULL)"

leading_eigens
all_sta = Float64[]; all_leading_collect = Vector{Vector{ComplexF64}}()
for j in 1: 6
        circ_leading_H = ComplexF64[];
        for i in 1:15
            push!(circ_leading_H, leading_eigens[j])
        end 
    push!(all_sta, sum(real.(circ_leading_H) .< 0)/length(path)); push!(all_leading_collect, circ_leading_H)
end 


path = glob("Eff_iters*", "../data/Eff_p0_old/")
all_circ = Vector{Vector{ComplexF64}}()
for j in 1: num_temps
    # if (j-1) % 5 == 0
        circ_leading_H = ComplexF64[]
        for i in 1:length(path)
            @load path[i] all_leading diag_dominance
            push!(circ_leading_H, all_leading[j])
        end 
    push!(all_circ, circ_leading_H)
end 
[sum(real.(all_circ[t]) .< 0) for t in 1:num_temps] ./ length(path)