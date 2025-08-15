
function growth_MiCRM!(dx, x, p, t, i)

    # maintenance
    dx[i] = -x[i]*p.m[i]

    # resoruce uptake and leakage
    for α = 1:p.M
        dx[i] += x[α + p.N] * x[i] * p.u[i, α]
        for β=1:p.M
            dx[i] += -x[i]*x[α + p.N]*p.u[i, α]*p.l[i, α, β]
        end
    end
end

function supply_MiCRM!(dx, x, p, t, α)
    if p.input_type == "constant"
        dx[α + p.N] = p.ρ[α]
    elseif p.input_type == "leaching"
        dx[α + p.N] = p.ρ[α] - (x[α + p.N] * p.ω[α])
    elseif p.input_type == "chemostat"
        dx[α + p.N] = p.ω[α] * (p.Kc[α] - x[α + p.N])
    elseif p.input_type == "self-renewing"
        r_self = fill(1,p.M) .* exp.((-0.32./0.0000862) * ((1/p.T)-(1/283.15)))./(1 .+ (0.32./(3.5 .- 0.32)) .* exp.(3.5/0.0000862 * (1 ./(273.15+25) .- 1/p.T)))
        dx[α + p.N] = (r_self[α] * x[α + p.N])/p.Kc[α] * (p.Kc[α] - x[α + p.N])
    end 
end

# x = sol.u[length(sol.t)-65]
# p.ρ[α]
# p.ρ[α] - (x[α + p.N] * p.ω[α])
# p.ω[α] * (p.Kc[α] - x[α + p.N])
# (r_self[α] * x[α + p.N])/p.Kc[α] * (p.Kc[α] - x[α + p.N])

function depletion_MiCRM!(dx, x, p, t, i, α)

    # uptake
    dx[α + p.N] += -(x[α + p.N] * p.u[i, α] * x[i])

    # leakage
    for β = 1:p.M
        dx[α + p.N] += x[β + p.N] * x[i] * p.u[i, β] * p.l[i, β, α]
    end
end

function dx!(dx, x, p, t;
    growth!::Function = growth_MiCRM!,
    supply!::Function = supply_MiCRM!,
    depletion!::Function = depletion_MiCRM!)

    for i = 1:p.N
        # reset derivatives
        dx[i] = 0.0

        # if consumer is extant
        if x[i] > 1e-5
            # update dx of ith consumer
            growth!(dx, x, p, t, i)
        end
    end

    for α = 1:p.M
        # reset derivatives
        dx[p.N + α] = 0.0

        #supply term
        supply!(dx, x, p, t, α)

        # loop over consumers
        for i = 1:p.N
            if x[i] > 1e-5
                depletion!(dx, x, p, t, i, α)
            end
        end
    end

end

function dxx!(dx, x, p, t)

    for i =1:p.N
        dx[i] = 0.0
        dx[i] = -p.m[i]*x[i]

        for α = 1:p.M
            dx[i] += x[i]*x[α + p.N]*p.u[i, α]
            for β=1:p.M
                dx[i] += -x[i]*x[α + p.N]*p.u[i, α]*p.l[i, α, β]
            end
        end
    end
    for α = 1:p.M
        dx[α + p.N] = 0.0
        if p.input_type == "constant"
            dx[α + p.N] = p.ρ[α]
        elseif p.input_type == "leaching"
            dx[α + p.N] = p.ρ[α] - (x[α + p.N] * p.ω[α])
        elseif p.input_type == "chemostat"
            dx[α + p.N] = p.ω[α] * (p.Kc[α] - x[α + p.N])
        elseif p.input_type == "self-renewing"
            r_self = fill(1,p.M) .* exp.((-0.32./0.0000862) * ((1/p.T)-(1/283.15)))./(1 .+ (0.32./(3.5 .- 0.32)) .* exp.(3.5/0.0000862 * (1 ./(273.15+25) .- 1/p.T)))
            dx[α + p.N] = (r_self[α] * x[α + p.N])/p.Kc[α] * (p.Kc[α] - x[α + p.N])
        end 


        for i=1:p.N
            dx[α + p.N] += -p.u[i, α]*x[α+p.N]*x[i]
            for β=1:p.M
                dx[α + p.N] += x[β + p.N] * x[i] * p.u[i, β] * p.l[i, β, α]
            end
        end
    end
end
