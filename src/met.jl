
abstract type AbstractMetParam end


@with_kw mutable struct Met{T<:AbstractMetParam}
    year::Int = 0
    obs::T = T()
    errors::Vector{Float64} = []
    values::Vector{Float64} = []
end


# One Met per year (row) of the stats table `df`, in year order. Field f of
# T is read from the columns "<f>_<kw>0" to "<f>_<kw>12" (whole year, then
# each month), or "<f>0" to "<f>12" when `kw` is empty.
function create_mets(::Type{T}, df::AbstractDataFrame, kw::AbstractString="") where T<:AbstractMetParam
    sfx = isempty(kw) ? "" : "_$(kw)"
    cols = [f => ["$(f)$(sfx)$(i)" for i ∈ 0:12] for f ∈ fieldnames(T)]
    map(eachrow(sort(df, :year))) do r
        obs = T(; (f => Float64[r[c] for c ∈ cs] for (f, cs) ∈ cols)...)
        Met{T}(year=r.year, obs=obs)
    end
end
