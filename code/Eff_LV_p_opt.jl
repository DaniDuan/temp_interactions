"""
Function for calculating effective Lotka-Volterra parameters.
It is necessary to supply the parameters p and an ODEProbelm
solution containing equilibrium values for the system in
question. verbose = true allows the function to return the
partial derivative vector as well as the A matrix.

"""

function Eff_LV_params(; p, sol) #, verbose = false

## Parameters are unpacked from dictionary and loaded into variables
    M, N, l, ρ, ω, m, u, λ = p.M, p.N, p.l, p.ρ, p.ω, p.m, p.u, p.λ

## We define a Kroenecker delta function
    δ(x, y) = ==(x, y)

## Equilibrium values are loaded into their respective vectors
    Ceq = sol[1:N, end]
    Req = sol[(N+1):(N+M), end]

## LV parameters are initialized
    A = zeros(M, M) # A matrix
    ∂R = zeros(M, N) # Partial derivative vector
    ℵ = zeros(N, N) # Interaction matrix
    O = zeros(N) # Dummy variable for calculating r
    P = zeros(N) # Dummy variable for calculating r
    r = zeros(N) # Effective growth rates

## Calculating the A matrix from MiCRM parameters and equilibrium solutions
    A = [(-ω[α] + sum(l[i, α, β] * u[i, β] * Ceq[i] - u[i, β] * Ceq[i] * δ(α, β) for i in 1:N)) for α in 1:M, β in 1:M]
    
## Precompute repeated values
    invA = inv(A)
    A_thing = u .* (1 .- λ) 

## Calculating partial derivatives of eq resources with respect to consumers
    ∂R = [sum(invA[α, β] * u[j, β] * Req[γ] * (δ(β, γ) - l[j, β, γ]) for β in 1:M, γ in 1:M) for α in 1:M, j in 1:N]

## Calculating interaction matrix
    ℵ = [sum(A_thing[i, α] * ∂R[α, j] for α in 1:M) for i in 1:N, j in 1:N]

## Calculating components of the intrinsic growth rates
    O = [sum(A_thing[i, α] * Req[α] for α in 1:M) for i in 1:N]

    P = [dot(ℵ[i, :], Ceq) for i in 1:N]
    
## Calculating effective intrinsic growth rates
    r = O .- P .- m

    ## Check verbose value and return corresponding parameter dictionary
    # if verbose == false
        return (ℵ = ℵ, r = r, N = N)
    # else
    #     return (ℵ = ℵ, r = r, N = N, ∂R = ∂R, A = A)
    # end
end


# @time p_lv = Eff_LV_params(p=p, sol=sol);
