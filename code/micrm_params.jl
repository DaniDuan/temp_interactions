

def_m(N, M, kw) = ones(N)
def_ρ(N, M, kw) = ones(M)
def_ω(N, M, kw) = ones(M)
def_u(N, M, kw) = copy(rand(Distributions.Dirichlet(M, 1.0), N)')

function def_l(N, M, L)
    Random.seed!(6)
    l = zeros(N, M, M)
    ϕ = fill(1.0, M)
    dD = Dirichlet(ϕ[:])
    for i = 1:N
        for α = 1:M
            l[i, α, :] = rand(dD) * L[i]
        end
    end
    return l
end

function generate_params(N, M; f_m=def_m, f_ρ=def_ρ, f_ω=def_ω, f_u=def_u, f_l=def_l, kwargs...)
    kw = Dict{Symbol, Any}(kwargs) # for temp, include: T, ρ_t, Tr, Ed, L

    tt, B, E, Tp = temp_trait(N, kw)
    push!(kw, :tt => tt)

    # consumers
     m = f_m(N, M, kw)
     u = f_u(N, M, kw)

     # Adding in the trade-off between growth and yield (0.1 to 0.7 for leakage & -2 to 4 for log(u_i))
     #  L = 0.1 .* log.(sum(u,dims = 2)) .+ 0.3
    # L = (log.(sum(u,dims = 2)) .+ 2) ./ 6 # -2 < min(log(u_i)), 4 > max(log(u_i)) 
     
     l = f_l(N, M, L)

    #  ## Calculating total leakage of consumers per resource
    # for i in 1:N
    #     for α in 1:M
    #         λ[i, α] = sum(l[i, α, :])
    #     end
    # end
    λ = reshape(sum(l , dims = 3), N, M)

     # resources
     ρ = f_ρ(N, M, kw)
     ω = f_ω(N, M, kw)

     kw_nt = (; kwargs...)
     p_nt = (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ, L = L, B, E, Tp)
    #  p_nt = (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ)

     out = Base.merge(p_nt, kw_nt)

     return out #(N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ)
 end

# p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
# # p.Eϵ
# u = p.u 
# u_sum = log.(sum(u,dims = 2))
# L = 0.1 .* u_sum .+ 0.3
# (u_sum .+2 ) ./ 6