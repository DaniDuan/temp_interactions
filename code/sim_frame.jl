# cd("WORKING_DIRECTORY")

# Load libraries
using Distributions
using LinearAlgebra
using DifferentialEquations
# using Plots, StatsPlots
using Sundials
using Parameters
using CSV, DataFrames
using CairoMakie
using LsqFit 

# Include simulation code files
include("micrm_params.jl") # Contains function gereate_params with default sampling scheme

include("dx_v2.jl") # Defines differential equations for MiCRM integration use dxx instead of dx

include("LV_dx.jl") # Defines LV differential equatons, use LV_dx

include("temp.jl")

include("EFF_LV_p_opt.jl")

####################################################################################################################################################################################
function F_m(N, M, kw)
    if haskey(kw, :T)
        m = kw[:tt][:,2]
    else 
        m = fill(0.2, N)
    end 
    return m
end

function F_ρ(N, M, kw)
    ρ = fill(1, M)
    return ρ
end

function F_ω(N, M, kw)
    ω = fill(0.0, M)
    return ω
end

function F_u(N, M, kw)
    if haskey(kw, :T)
        u_sum = kw[:tt][:,1]
    else 
        u_sum = fill(2.5, M) 
    end
    diri = transpose(rand(Dirichlet(ones(M)),N))
    u = diri.*u_sum
    return u
end
