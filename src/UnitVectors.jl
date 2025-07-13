module UnitVectors

using Random
using StaticArrays
using LinearAlgebra

export UnitVector

struct UnitVector <: FieldVector{3, Float64}
    x::Float64
    y::Float64
    z::Float64

    UnitVector() = new(0, 0, 1)
    UnitVector(x, y, z) =
        norm((x, y, z)) ≈ 1.0 ? new(x, y, z) : error("Non-unit unit vector")
    UnitVector((x, y, z)::Tuple) =
        norm((x, y, z)) ≈ 1.0 ? new(x, y, z) : error("Non-unit unit vector")
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{UnitVector})
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return UnitVector(cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end

end