
# How each temperature and wind statistic is measured from daily values
const STAT_FNS = (mean=mean, sd=x -> max(0.01, std(x)), rlag=acf1, skew=skewness)


# Statistics of one year of daily `data`: whole year, then each month
function collect_met(::Type{T}, year::Int, data) where T<:Union{Temp,Wind}
    parts = year_and_months(data, year)
    T(; (f => getproperty(STAT_FNS, f).(parts) for f ∈ fieldnames(T))...)
end


function collect_met(::Type{Rain}, year::Int, data)
    totrain, pww, pwd = partition_rain(year, data)
    Rain(totrain=totrain, pww=pww, pwd=pwd)
end


# Variables of the stats file, in order: column of daily weather, its
# statistics, and the suffix of their column names
const STATS_VARS = (("tmin", Temp, "_tmin"), ("tmax", Temp, "_tmax"),
                    ("wind", Wind, "_wind"), ("rain", Rain, ""))


"""
    weather_stats(df)

Monthly and annual statistics of the daily weather `df` (columns `year` and
any of `tmin`, `tmax`, `wind`, `rain`; whole years only), one row per year,
in the layout of the stats file.
"""
function weather_stats(df::AbstractDataFrame)
    vars = [v for v ∈ STATS_VARS if v[1] ∈ names(df)]
    rows = map(collect(groupby(df, :year))) do g
        year = g.year[1]
        cols = Pair{Symbol,Any}[:year => year]
        for (col, T, sfx) ∈ vars
            met = collect_met(T, year, Vector{Float64}(g[!, col]))
            for f ∈ fieldnames(T), (i, v) ∈ enumerate(getproperty(met, f))
                push!(cols, Symbol("$(f)$(sfx)$(i-1)") => v)
            end
        end
        (; cols...)
    end
    DataFrame(rows)
end


# Writes the stats file: latitude line, header, then one row per year
function write_stats(fname::AbstractString, lat, stats::AbstractDataFrame)
    open(fname, "w") do fout
        println(fout, lat)
        println(fout, join(names(stats), ","))
        for r ∈ eachrow(stats)
            println(fout, join(r, ","))
        end
    end
    fname
end


function create_data_file(wthrfname::AbstractString, inputfname::AbstractString)
    wthr = csv2df(wthrfname)
    fullpath_inputfname = joinpath(dirname(wthrfname), inputfname)
    write_stats(fullpath_inputfname, wthr.lat, weather_stats(wthr.df))
end
