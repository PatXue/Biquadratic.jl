const SpinVector = SVector{3, Float64}
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

# Sum spins of position (x, y)'s nearest neighbors
nn_sum(A::AbstractArray, x, y) = A[x+1, y] + A[x, y+1] + A[x-1, y] + A[x, y-1]
# Sum spins of (x, y)'s next nearest NE-SW neighbors
nnna_sum(A::AbstractArray, x, y) = A[x+1, y+1] + A[x-1, y-1]
# Sum spins of (x, y)'s next nearest NW-SE neighbors
nnnb_sum(A::AbstractArray, x, y) = A[x+1, y-1] + A[x-1, y+1]
