
abstract type AbstractMetParam end


@with_kw mutable struct Met{T<:AbstractMetParam}
    year::Int = 0
    obs::T = T()
    est::T = T()
    errors::Vector{Float64} = []
    values::Vector{Float64} = []
end
