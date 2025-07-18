module Biquadratic

import Random.AbstractRNG
import Random.default_rng

using Carlo
using HDF5
using LinearAlgebra
using PeriodicArrays
using Random
using StaticArrays

const SpinVector = SVector{3, Float64}
# Note: Using temperature in units of energy (k_B = 1)
# All energy units are in terms of J2a (best to set J2a = 1)
struct MC{AlgType} <: AbstractMC
    T::Float64   # Temperature
    J1::Float64  # Nearest neighbor coupling energy
    J2a::Float64 # Next-nearest neighbor coupling energy (NE-SW direction)
    J2b::Float64 # Next-nearest neighbor coupling energy (NW-SE direction)
    K::Float64   # Biquadratic coupling energy

    spins::PeriodicMatrix{SpinVector}
end

function MC{AlgType}(T=0.5, J1=0.1, J2a=1.0, J2b=-1.0, K=0.1,
                 Lx::Int=20, Ly::Int=20) where {AlgType}
    MC{AlgType}(T, J1, J2a, J2b, K, fill(zeros(SpinVector), (Lx, Ly)))
end

function MC{AlgType}(params::AbstractDict) where {AlgType}
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    J1 = params[:J1]
    J2a = params[:J2a]
    J2b = params[:J2b]
    K = params[:K]
    return MC{AlgType}(T, J1, J2a, J2b, K, Lx, Ly)
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{SpinVector})
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return SpinVector(cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end

function init_eag!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, y = Tuple(I)
        spin_sign = mod(x+y, 4) < 2 ? 1.0 : -1.0
        spins[I] = spin_sign * SVector(1.0, 0.0, 0.0)
    end
end

function init_orth!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, y = Tuple(I)
        θ = π/2 * (x + y)
        spins[I] = SVector(cos(θ), sin(θ), 0.0)
    end
end

function Carlo.init!(mc::MC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = params[:init_type]
    if init_type == :const
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :orth
        init_orth!(mc.spins)
    elseif init_type == :eag
        init_eag!(mc.spins)
    else
        rand!(ctx.rng, mc.spins)
    end
    return nothing
end

# Sum spins of position (x, y)'s nearest neighbors
nn_sum(A::AbstractArray, x, y) = A[x+1, y] + A[x, y+1] + A[x-1, y] + A[x, y-1]
# Sum spins of (x, y)'s next nearest NE-SW neighbors
nnna_sum(A::AbstractArray, x, y) = A[x+1, y+1] + A[x-1, y-1]
# Sum spins of (x, y)'s next nearest NW-SE neighbors
nnnb_sum(A::AbstractArray, x, y) = A[x+1, y-1] + A[x-1, y+1]

# Calculate the energy at a lattice site (x, y) if it had spin s
function energy(mc::MC{:Metropolis}, s::SpinVector, x, y)
    nn = nn_sum(mc.spins, x, y)
    nnna = nnna_sum(mc.spins, x, y)
    nnnb = nnnb_sum(mc.spins, x, y)
    H0 = s ⋅ (mc.J1 * nn + mc.J2a * nnna + mc.J2b * nnnb)
    biquad = mc.K * ((s ⋅ mc.spins[x-1, y])^2 + (s ⋅ mc.spins[x+1, y])^2)
    return H0 + biquad
end

function Carlo.sweep!(mc::MC{:Metropolis}, rng::AbstractRNG=default_rng())
    Lx, Ly = size(mc.spins)
    for _ in 1:length(mc.spins)
        # Select site for spin change
        x = rand(rng, 1:Lx)
        y = rand(rng, 1:Ly)
        old_s = mc.spins[x, y]
        # Propose new spin vector
        new_s = rand(rng, SpinVector)
        ΔE = energy(mc, new_s, x, y) - energy(mc, old_s, x, y)

        # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
        prob = exp(-ΔE / mc.T)
        if prob >= 1.0 || rand(rng) < prob
            mc.spins[x, y] = new_s
        end
    end
    return nothing
end

function Carlo.sweep!(mc::MC, ctx::Carlo.MCContext)
    Carlo.sweep!(mc, ctx.rng)
end

# Calculate the energy contribution of a site (x, y), considering only half of
# its bonds (avoids double counting when calculating total energy)
function half_energy(mc::MC{:Metropolis}, x, y)
    s = mc.spins[x, y]
    nn = mc.spins[x+1, y] + mc.spins[x, y+1]
    nnna = mc.spins[x+1, y+1]
    nnnb = mc.spins[x+1, y-1]
    H0 = s ⋅ (mc.J1 * nn + mc.J2a * nnna + mc.J2b * nnnb)
    biquad = mc.K * (s ⋅ mc.spins[x+1, y])^2
    return H0 + biquad
end

function Carlo.measure!(mc::MC{:Metropolis}, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    # Magnetization per lattice site
    mag = norm(sum(mc.spins)) / N
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Mag4, mag^4)

    # Energy per lattice site
    energy = 0.0
    # Averaged adjacent dot products
    Dx0 = Dy0 = 0.0
    Dxπ = Dyπ = 0.0
    # Spin current
    P = zeros(3)
    x_hat = [1.0, 0.0, 0.0]
    for y in 1:Ly
        for x in 1:Lx
            energy += half_energy(mc, x, y)

            s = mc.spins[x, y]
            sx = mc.spins[x+1, y]
            sy = mc.spins[x, y+1]
            x_dot = s ⋅ sx
            y_dot = s ⋅ sy
            Dx0 += x_dot
            Dy0 += y_dot
            Dxπ += x_dot * (-1)^(x+y)
            Dyπ += y_dot * (-1)^(x+y)

            P += x_hat × (s × sx)
        end
    end
    energy /= N
    measure!(ctx, :Energy, energy)
    measure!(ctx, :Energy2, energy^2)
    Dx0 /= N
    Dy0 /= N
    measure!(ctx, :Dx0, Dx0)
    measure!(ctx, :Dy0, Dy0)
    Dxπ /= N
    Dyπ /= N
    measure!(ctx, :Dxπ, abs(Dxπ))
    measure!(ctx, :Dyπ, abs(Dyπ))
    P /= N
    measure!(ctx, :P, norm(P))

    return nothing
end

function Carlo.register_evaluables(::Type{MC{AlgType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType}
    T = params[:T]
    J2a = params[:J2a]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N * J2a/T * (mag2 - mag^2)
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
    raw_spins = read(in, "spins")
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_spins)
    return nothing
end

end