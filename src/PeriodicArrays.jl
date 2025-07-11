module PeriodicArrays

export PeriodicArray, PeriodicMatrix

struct PeriodicArray{T, N} <: AbstractArray{T, N}
    array::Array{T, N}
end
const PeriodicMatrix{T} = PeriodicArray{T, 2}

Base.size(A::PeriodicArray) = size(A.array)
Base.convert(::Type{PeriodicArray{T, N}}, A::AbstractArray{T, N}) where {T, N} =
    PeriodicArray(A)

function Base.getindex(A::PeriodicArray, I::Vararg{Int, N}) where {N}
    indices = map(mod1, I, size(A.array))
    return getindex(A.array, indices...)
end
function Base.setindex!(A::PeriodicArray, v, I::Vararg{Int, N}) where {N}
    indices = map(mod1, I, size(A.array))
    setindex!(A.array, v, indices...)
end

end