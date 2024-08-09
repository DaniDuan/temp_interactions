include("./sim_frame.jl");

N=100
M=50
L = 0.3
### Temp params 
num_temps = 38
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


all_ℵii = Vector{Vector{Float64}}(); all_ℵij = Vector{Vector{Float64}}(); all_up_ℵij = Vector{Vector{Float64}}(); all_low_ℵij = Vector{Vector{Float64}}()
for i in range(0, stop = num_temps-1, length = num_temps)
    T = 273.15 + i
    Random.seed!(6)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    # println(bm .> 1.0e-7)
    sur = (1:N)[bm .> 1.0e-7]
    N_s = length(sur)
    p_lv = Eff_LV_params(p=p, sol=sol);
    # number of species with r>0 at equilibium 
    sur_ℵ = p_lv.ℵ
    ℵii = diag(sur_ℵ)
    # ℵij = vec(sum(reshape(sur_ℵ,Int(sqrt(length(sur_ℵ))),Int(sqrt(length(sur_ℵ)))), dims = 2).- diag(reshape(sur_ℵ,Int(sqrt(length(sur_ℵ))),Int(sqrt(length(sur_ℵ))))))
    ℵ = reshape(sur_ℵ,Int(sqrt(length(sur_ℵ))),Int(sqrt(length(sur_ℵ))))
    ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i != j]
    up_ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i > j]
    low_ℵij = [ℵ[i,j] for i in 1:N for j in 1:N if i < j]
    push!(all_ℵii, ℵii); push!(all_ℵij, ℵij); push!(all_up_ℵij, up_ℵij); push!(all_low_ℵij, low_ℵij)
    print(i, " °C Complete, ", "α ",mean(ℵii),"\n") 

end 

Temp_rich = range(0, num_temps-1, length = num_temps)

using LsqFit, GLM
k = 0.0000862 # Boltzman constant
temp_SS(T, params) = params[1] .* exp.((-params[2]./k) * ((1 ./T) .-(1/Tr)))./(1 .+ (params[2]./(params[4] .- params[2])) .* exp.(params[4]/k * (1 ./params[3] .- 1 ./T)))

Random.seed!(6); 
x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
temp = collect(Temp_rich .+273.15)

f1 = Figure(resolution = (1200, 1200));
f2 = Figure(resolution = (1200, 1200));
# f3 = Figure(resolution = (1200, 1200));
# f4 = Figure(resolution = (1200, 1200));
ax1 = Axis(f1[1,1], xlabel = "Temperature (°C)", ylabel = "αii", xlabelsize = 50, ylabelsize = 50)
ax2 = Axis(f2[1,1], xlabel = "Temperature (°C)", ylabel = "αij", xlabelsize = 50, ylabelsize = 50)

Temp_α = zeros(Float64, N, 11)
for i in 1:N
    αii = [all_ℵii[t][i] for t in 1:num_temps]
    Tpii = Temp_rich[argmax(abs.(αii))] + 273.15
    Eii = maximum(diff(log.(abs.(αii)))./diff(x))
    init_ii = [abs(αii[Int(Tr-273.15+1)]), Eii, Tpii, 3.5]

    αij = [all_ℵij[t][i] for t in 1:num_temps]
    Tpij = Temp_rich[argmax(abs.(αij))] + 273.15
    Eij = maximum(diff(log.(abs.(αij)))./diff(x))
    init_ij = [abs(αij[Int(Tr-273.15+1)]), Eij, Tpij, 3.5]

    fit_ii = curve_fit(temp_SS, temp, abs.(αii), init_ii)
    fit_ij = curve_fit(temp_SS, temp, abs.(αij), init_ij)
    # ax1 = Axis(f1[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    # ax2 = Axis(f2[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    # ax3 = Axis(f3[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    # ax4 = Axis(f4[Int(floor((i-1)/10+1)),Int((i-1) % 10+1)], ygridvisible = false, xgridvisible = false)
    # lines!(ax1, Temp_rich, abs.(αii), color = ("#285C93", 1), linewidth = 1)
    # lines!(ax2, Temp_rich, abs.(temp_SS(temp, fit_ii.param)), color = ("#E17542", 1), linewidth = 1)
    # lines!(ax3, Temp_rich, abs.(αij), color = ("#285C93", 1), linewidth = 1)
    # lines!(ax4, Temp_rich, abs.(temp_SS(temp, fit_ij.param)), color = ("#E17542", 1), linewidth = 1)
    scatter!(ax1, Temp_rich, log.(abs.(αii)), color = "#285C93", alpha = 0.5)
    scatter!(ax2, Temp_rich, log.(abs.(αij)), color = "#E17542", alpha = 0.5)

    Temp_α[Int(i),:] = vcat(Int(i),fit_ii.param, AIC(fit_ii, num_temps), fit_ij.param, AIC(fit_ij, num_temps))
end 
col_names_α = ["species","Bii", "Eii", "Tpii", "Edii", "AIC_ii", "Bij", "Eij", "Tpij", "Edij", "AIC_ij"];
Temp_α = DataFrame(Temp_α, col_names_α);
f1
# CSV.write("../data/Tα.csv", Temp_α, writeheader=true)

temp = collect(Temp_rich .+273.15)
temp_all = repeat(temp ,inner = 100)
allii = vcat(all_ℵii...)


ii = [mean(log.(abs.(all_ℵii[t][i])) for i in 1:N) for t in 1:num_temps]
Tpii_m = Temp_rich[argmax(ii)]
Bii_m = ii[Int(Tr-273.15+1)]
Eii_up = maximum(diff(ii)./diff(x))
Eiir = range(0, Eii_up, 50000)

Eii_all = zeros(Float64, 50000, 4)
for i in 1: 50000
    Bii = exp(rand(Normal(Bii_m, 5)))
    Eii = rand(Uniform(0, Eii_up))
    Tpii = 273.15 .+ rand(Normal(Tpii_m, 10))
    init_ii = [Bii, Eii, Tpii, 3.5]
    fit_ii = curve_fit(temp_SS, temp_all, abs.(allii), init_ii)
    AICii = AIC(fit_ii, num_temps * N)
    Eii_all[Int(i),:] = [Bii, Eii, Tpii, AICii]
    if i%500 == 0
        print(i/50000, "\n")
    end 
end 

Bii = Eii_all[:,1][argmin(Eii_all[:,4])]
Eii = Eii_all[:,2][argmin(Eii_all[:,4])]
Tpii = Eii_all[:,3][argmin(Eii_all[:,4])]

init_ii = [Bii, Eii, Tpii, 3.5]

fit_ii = curve_fit(temp_SS, temp_all, abs.(allii), init_ii)
AIC(fit_ii, num_temps * N)
fit_ii.param

f1 = Figure(resolution = (1200, 1200));
ax1 = Axis(f1[1,1], xlabel = "Temperature (°C)", ylabel = "αii", xlabelsize = 50, ylabelsize = 50)
scatter!(ax1, temp_all .- 273.15, log.(abs.(allii)), color = "#285C93", alpha = 0.5)
lines!(ax1, Temp_rich, log.(abs.(temp_SS(temp, fit_ii.param))), color = ("#E17542", 1), linewidth = 1)
lines!(ax1, Temp_rich, ii, color = ("#285C93", 1), linewidth = 1)
f1