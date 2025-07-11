module PeriodicArrays

struct PeriodicArray{T, N} <: AbstractArray{T, N}
    A::Array{T, N}
end
const PeriodicMatrix{T} = PeriodicArray{T, 2}

end