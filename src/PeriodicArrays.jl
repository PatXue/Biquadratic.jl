module PeriodicArrays

export PeriodicArray, PeriodicMatrix

struct PeriodicArray{T, N} <: AbstractArray{T, N}
    array::Array{T, N}
end
const PeriodicMatrix{T} = PeriodicArray{T, 2}

Base.size(A::PeriodicArray) = size(A.array)
function Base.getindex(A::PeriodicArray{T, N}, I::Vararg{Int, N}) where {T, N}
    indices = map(mod1, I, size(A.array))
    return getindex(A.array, indices...)
end

end