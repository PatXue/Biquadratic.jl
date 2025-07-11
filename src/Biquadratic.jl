module Biquadratic

include("PeriodicArrays.jl")

import Random.AbstractRNG
import Random.default_rng

using .PeriodicArrays
using Carlo
using HDF5
using LinearAlgebra
using StaticArrays

# Note: Using temperature in units of energy (k_B = 1)
struct MC <: AbstractMC
    T::Float64   # Temperature
    J1::Float64  # Nearest neighbor coupling energy
    J2a::Float64 # Next-nearest neighbor coupling energy (NE-SW direction)
    J2b::Float64 # Next-nearest neighbor coupling energy (NW-SE direction)
    K::Float64   # Biquadratic coupling energy

    spins::PeriodicMatrix{SVector{3, Float64}}
end

function MC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    J1 = params[:J1]
    J2a = params[:J2a]
    J2b = params[:J2b]
    K = params[:K]
    return MC(T, J1, J2a, J2b, K, fill(SVector(0, 0, 0), (Lx, Ly)))
end

"""
    rand_vector(vec, [rng = default_rng()])

Fill the first 3 elements of vec with a 3D unit vector, generated uniformly
"""
function rand_vector!(vec, rng::AbstractRNG=default_rng())
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)

    vec[1] = cos(ϕ)sin(θ)
    vec[2] = sin(ϕ)sin(θ)
    vec[3] = cos(θ)
    return nothing
end

function Carlo.init!(mc::MC, ctx::Carlo.MCContext, params::AbstractDict)
    vec::Vector{Float64} = [0, 0, 1]
    for slice in eachslice(mc.spins, dims=(1, 2))
        if params[:rand_init]
            rand_vector!(vec, ctx.rng)
        end
        slice .= vec
    end

    return nothing
end

function sum_adj(M, (x, y))
    Lx, Ly = size(M, 1), size(M, 2)
    return M[mod1(x-1, Lx), y, :] + M[x, mod1(y-1, Ly), :] +
           M[mod1(x+1, Lx), y, :] + M[x, mod1(y+1, Ly), :]
end

function Carlo.sweep!(mc::MC, rng::AbstractRNG=default_rng())
    Lx, Ly = size(mc.spins, 1), size(mc.spins, 2)
    for _ in 1:length(mc.spins)
        # Select site for spin change
        x = rand(rng, 1:Lx)
        y = rand(rng, 1:Ly)

        # Sum of nearest neighbors' spins
        H = mc.J * sum_adj(mc.spins, (x, y)) / mc.T
        unit_H = H / norm(H)
        H⊥ = nullspace(unit_H')

        # Randomly generate new θ and ϕ according to Boltzmann distribution
        # (relative to H)
        cosθ = log1p(rand(rng) * (exp(2norm(H)) - 1)) / norm(H) - 1
        ϕ = 2π * rand(rng)

        ϕ_comp = sqrt(1 - cosθ^2) * (cos(ϕ)H⊥[:, 1] + sin(ϕ)H⊥[:, 2])
        mc.spins[x, y, :] .= cosθ * unit_H + ϕ_comp
    end
    return nothing
end

function Carlo.sweep!(mc::MC, ctx::Carlo.MCContext)
    Carlo.sweep!(mc, ctx.rng)
end

function Carlo.measure!(mc::MC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins, 1), size(mc.spins, 2)
    N = Lx * Ly
    # Magnetization per lattice site
    mag = norm(sum(mc.spins, dims=(1, 2))) / N
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Mag4, mag^4)

    # Energy per lattice site
    energy = 0.0
    for x in 1:Lx
        for y in 1:Ly
            energy += -mc.J * mc.spins[x, y, :] ⋅
                (mc.spins[mod1(x+1, Lx), y, :] + mc.spins[x, mod1(y+1, Ly), :])
        end
    end
    energy /= N
    energy += -mc.H * mag
    measure!(ctx, :Energy, energy)
    measure!(ctx, :Energy2, energy^2)

    return nothing
end

function Carlo.register_evaluables(
    ::Type{MC}, eval::AbstractEvaluator, params::AbstractDict
)
    T = params[:T]
    J = params[:J]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N * J/T * (mag2 - mag^2)
    end

    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end

    return nothing
end

function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["spins"] = mc.spins
    return nothing
end
function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.spins .= read(in, "spins")
    return nothing
end
using Carlo
using HDF5

end
