
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
        # E_r_self = 0.8; Tr = 283.15; k = 0.0000862

        #### For self-renewing resource with temperature dependent growth 
        # E_r_self = 0.8; Tpk_r_self = 273.15+25; Ed = 3.5; Tr = 283.15; k = 0.0000862
        # r_self = fill(1,p.M) * exp((-E_r_self/k) * ((1/p.T)-(1/Tr)))/(1 + (E_r_self/(Ed - E_r_self)) * exp(Ed/k * (1 /Tpk_r_self - 1/p.T)))
        # Kc_self = p.Kc * exp((0.8/k) * ((1/p.T)-(1/Tr)))
        # dx[α + p.N] = (r_self[α] * x[α + p.N])/Kc_self[α] * (Kc_self[α] - x[α + p.N])
        
        dx[α + p.N] = (fill(1,p.M)[α] * x[α + p.N])/p.Kc[α] * (p.Kc[α] - x[α + p.N])
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
            # E_r_self = 0.8; Tr = 283.15; k = 0.0000862
            # r_self = fill(1,p.M) * exp((-E_r_self/k) * ((1/p.T)-(1/Tr)))
            E_r_self = 0.8; Tpk_r_self = 273.15+25; Ed = 3.5; Tr = 283.15; k = 0.0000862
            r_self = fill(1,p.M) * exp((-E_r_self/k) * ((1/p.T)-(1/Tr)))/(1 + (E_r_self/(Ed - E_r_self)) * exp(Ed/k * (1 /Tpk_r_self - 1/p.T)))
            Kc_self = p.Kc * exp((0.8/k) * ((1/p.T)-(1/Tr)))
            dx[α + p.N] = (r_self[α] * x[α + p.N])/Kc_self[α] * (Kc_self[α] - x[α + p.N])
            # dx[α + p.N] = (r_self[α] * x[α + p.N])/p.Kc[α] * (p.Kc[α] - x[α + p.N])
        end 

        for i=1:p.N
            dx[α + p.N] += -p.u[i, α]*x[α+p.N]*x[i]
            for β=1:p.M
                dx[α + p.N] += x[β + p.N] * x[i] * p.u[i, β] * p.l[i, β, α]
            end
        end
    end
end
