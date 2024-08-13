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
    B_m = α[Int(Tr-273.15+1)]
    clean = α[.!ismissing.(α)]
    x_c = x[.!ismissing.(α)]
    E_up = maximum(diff(clean)./diff(x_c))
    Ed = 3.5
    temp = collect(Temp_rich .+273.15)

    allα = abs.(vcat(αs...))
    nm_index = .!ismissing.(allα)
    temp_all = vcat([repeat([temp[t]], length(αs[t])) for t in 1:nt]...)
    temp_all = temp_all[nm_index]
    allα = allα[nm_index]
    return Nα, B_m, E_up, T_m, Ed,temp_all, allα
end 
using LsqFit

using LsqFit

function calculate_r2(model, xdata, ydata)
    y_pred = temp_SS(xdata, model.param)
    ss_res = sum((ydata .- y_pred).^2)
    ss_tot = sum((ydata .- mean(ydata)).^2)
    r_square = 1 - ss_res / ss_tot
    return r_square
end

function try_params(αs, nt, iters = 1000)
    Nα, B_m, E_up, T_m, Ed_m, temp_all, allα = get_init_param(αs, nt)
    collection = zeros(Float64, iters, 5)
    Ed = Ed_m
    Tp = T_m +273.15
    for i in 1:iters
        B = exp(rand(Normal(B_m, 5)))
        E = rand(Uniform(0, E_up))
        # rand_0 = truncated(Normal(0,0.5), 0, Inf)
        # rand35 = truncated(Normal(3.5,1), 0, Inf)
        # Ed = rand([rand(rand_0), rand(rand35)])
        Ed = 3.5
        init = [B, E, Tp, Ed]
        try
            fit_ss = curve_fit(temp_SS, temp_all, allα, init)
            AIC_ss = AIC(fit_ss, Nα)
            collection[Int(i),:] = [B, E, Tp, Ed, AIC_ss]
        catch e 
            collection[Int(i),:] = [0.0, 0.0, 0.0, Ed, Inf]
        end
    end 
    B_in, E_in, T_in, Ed_in, AIC_in = collection[argmin(collection[:,5]),:]
    init_in = [B_in, E_in, T_in, Ed_in]
    return Nα, init_in, AIC_in, temp_all, allα
end 

