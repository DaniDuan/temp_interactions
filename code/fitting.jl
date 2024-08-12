include("./sim_frame.jl");

########### Defining functions #################
using LsqFit, GLM
temp_SS(T, params) = params[1] .* exp.((-params[2]./k) * ((1 ./T) .-(1/Tr)))./(1 .+ (params[2]./(params[4] .- params[2])) .* exp.(params[4]/k * (1 ./params[3] .- 1 ./T)))

countnonmiss(vec) = count(x -> !ismissing(x), vec)

function get_init_param(αs, nt)
    k = 0.0000862; Tr = 273.15+10
    Temp_rich = range(0, nt-1, length = nt)
    x = -1/k .* (1 ./(Temp_rich .+273.15) .- 1/Tr)
    α = [mean(log.(abs.(αs[t]))) for t in 1:nt]
    Nα = sum(countnonmiss.(αs))
    T_m = Temp_rich[argmax(α)]
    B_m = (α[Int(Tr-273.15)]+ α[Int(Tr-273.15+1)] + α[Int(Tr-273.15+2)])/3
    clean = α[.!ismissing.(α)]
    x_c = x[.!ismissing.(α)]
    E_up = maximum(diff(clean)./diff(x_c))
    Ed = -(clean[length(clean)] - clean[Int(T_m)])/(x_c[length(clean)] - x_c[Int(T_m)])
    if T_m < (nt-2) && Ed > 0
        Ed = Ed
    else 
        Ed = 3.5
    end 
    # E_mean = mean(diff(clean)./diff(x_c))
    temp = collect(Temp_rich .+273.15)

    allα = abs.(vcat(αs...))
    nm_index = .!ismissing.(allα)
    temp_all = vcat([repeat([temp[t]], length(αs[t])) for t in 1:nt]...)
    temp_all = temp_all[nm_index]
    allα = allα[nm_index]
    return Nα, B_m, E_up, T_m, Ed,temp_all, allα
end 

function try_params(αs, nt, iters = 1000)
    Nα, B_m, E_up, T_m, Ed_m, temp_all, allα = get_init_param(αs, nt)
    collection = zeros(Float64, iters, 5)
    for i in 1:iters
        B = exp(rand(Normal(B_m, 5)))
        E = rand(Uniform(0, E_up))
        Tp = T_m * rand() +273.15
        Ed = rand(Normal(Ed_m, Ed_m/2))
        init = [B, E, Tp, Ed]
        try
            fit_ss = curve_fit(temp_SS, temp_all, allα, init)
            AIC_ss = AIC(fit_ss, Nα)
            collection[Int(i),:] = [B, E, Tp, Ed, AIC_ss]
        catch e 
            collection[Int(i),:] = [0.0, 0.0, 0.0, 0.0]
        end
    end 
    B_in, E_in, T_in, Ed_in, AIC_in = collection[argmin(collection[:,5]),:]
    init_in = [B_in, E_in, T_in, Ed_in]
    return Nα, init_in, AIC_in, temp_all, allα
end 

