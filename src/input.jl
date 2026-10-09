
# How each temperature and wind statistic is measured from daily values
const STAT_FNS = (mean=mean, sd=x -> max(0.01, std(x)), rlag=acf1, skew=skew0)


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


# Dates of the rows of the daily weather `df`, after checking that they
# make whole years with every day once. Errors name the first problem.
function weather_dates(df::AbstractDataFrame)
    for c ∈ ("year", "month", "day")
        c ∈ names(df) || error("the daily weather needs a column `$(c)`")
    end
    dates = map(enumerate(zip(df.year, df.month, df.day))) do (i, (y, m, d))
        dt = try Date(Int(y), Int(m), Int(d)) catch; nothing end
        isnothing(dt) && error("data row $(i): $(y)-$(m)-$(d) is not a date")
        dt
    end
    seen = Set{Date}()
    for dt ∈ dates
        dt ∈ seen && error("$(dt) appears more than once in the daily weather")
        push!(seen, dt)
    end
    for y ∈ unique(Dates.year.(dates))
        n = count(dt -> Dates.year(dt) == y, dates)
        if n != daysinyear(y)
            gap = first(dt for dt ∈ Date(y):Day(1):Date(y, 12, 31) if dt ∉ seen)
            error("year $(y) has $(n) of its $(daysinyear(y)) days (the first missing " *
                  "is $(gap)); only whole years can be used")
        end
    end
    dates
end


# Errors if a weather column has a missing or non-numeric value, naming
# the first such date
function check_values(df::AbstractDataFrame, col::AbstractString, dates)
    i = findfirst(v -> !(v isa Real) || !isfinite(v), df[!, col])
    isnothing(i) || error("`$(col)` has no valid value on $(dates[i])")
end


"""
    weather_stats(df)

Monthly and annual statistics of the daily weather `df` (columns `year`,
`month`, `day` and any of `tmin`, `tmax`, `wind`, `rain`; whole years only,
in any row order), one row per year, in the layout of the stats file.
"""
function weather_stats(df::AbstractDataFrame)
    dates = weather_dates(df)
    order = sortperm(dates)
    df, dates = df[order, :], dates[order]
    vars = [v for v ∈ STATS_VARS if v[1] ∈ names(df)]
    foreach(v -> check_values(df, v[1], dates), vars)
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


# Allowed monthly values of each statistic in the stats file: a test, and
# what it requires
const STAT_RULES = Dict(
    "mean_t" => (isfinite, "a number"),
    "sd_t" => (>(0), "above 0"),
    "rlag_t" => (v -> abs(v) < 1, "between -1 and 1"),
    "skew_t" => (isfinite, "a number"),
    "mean_wind" => (>(0), "above 0"),
    "sd_wind" => (>(0), "above 0"),
    "rlag_wind" => (v -> abs(v) < 1, "between -1 and 1"),
    "totrain" => (>=(0), "0 or more"),
    "pww" => (v -> 0 <= v <= 1, "between 0 and 1"),
    "pwd" => (v -> 0 <= v <= 1, "between 0 and 1"))


# Checks the stats table `df` before generating from it: one row per year,
# every monthly column of each variable present, and every monthly value
# valid (STAT_RULES). Errors list the first problems found.
function check_stats(df::AbstractDataFrame)
    "year" ∈ names(df) || error("the stats file needs a column `year`")
    allunique(df.year) || error("the stats file has a year more than once")
    problems = String[]
    for (col, T, sfx) ∈ STATS_VARS
        any(occursin.(col, names(df))) || continue
        for f ∈ fieldnames(T), i ∈ 1:12
            name = "$(f)$(sfx)$(i)"
            if name ∉ names(df)
                push!(problems, "column $(name) is missing")
                continue
            end
            key = startswith(sfx, "_t") ? "$(f)_t" : "$(f)$(sfx)"
            ok, needs = STAT_RULES[key]
            for (y, v) ∈ zip(df.year, df[!, name])
                (v isa Real && isfinite(v) && ok(v)) ||
                    push!(problems, "$(name) of $(y) is $(v); it must be $(needs)")
            end
        end
    end
    isempty(problems) && return nothing
    shown = first(problems, 10)
    more = length(problems) > 10 ? "\n  … and $(length(problems) - 10) more" : ""
    error("the stats file has $(length(problems)) problem(s):\n  " * join(shown, "\n  ") * more)
end
