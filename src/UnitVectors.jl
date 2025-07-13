module UnitVectors

using Random
using LinearAlgebra

export UnitVector

struct UnitVector <: AbstractArray{Float64, 1}
    x::Float64
    y::Float64
    z::Float64

    UnitVector() = new(0, 0, 1)
    UnitVector(x, y, z) =
        norm((x, y, z)) ≈ 1.0 ? new(x, y, z) : error("Non-unit unit vector")
    UnitVector((x, y, z)::Tuple) =
        norm((x, y, z)) ≈ 1.0 ? new(x, y, z) : error("Non-unit unit vector")
    function UnitVector(v::AbstractVector)
        if length(v) != 3
            error("UnitVector must be 3D")
        end
        return UnitVector(v[1], v[2], v[3])
    end
end

Base.size(::UnitVector) = (3,)

const symtable::NTuple{3, Symbol} = (:x, :y, :z)
# Indexing functions
function Base.getindex(v::UnitVector, i::Int)
    if i < 1 || 3 < i
        Base.throw_boundserror(v, i)
    else
        return getfield(v, symtable[i])
    end
end
Base.IndexStyle(::Type{UnitVector}) = IndexLinear()

# Random generation
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{UnitVector})
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return UnitVector(cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end

end