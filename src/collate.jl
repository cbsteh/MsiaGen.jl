
function param_fields(obj::AbstractMetParam)
    est_fields = [isa(getproperty(obj, field), AbstractVector) ?
                 ["$(field)$(i)" for i ∈ 0:12] :
                 ["$(field)"] for field ∈ propertynames(obj)]
    collect(Iterators.flatten(est_fields))
end


function param_values(obj::AbstractMetParam)
    est_values = [isa(getproperty(obj, field), AbstractVector) ?
                 [getproperty(obj, field)[i] for i ∈ 1:13] :
                 [getproperty(obj, field)] for field ∈ propertynames(obj)]
    collect(Iterators.flatten(est_values))
end


function collect_mettype(lst, prefix, mettype)
    df = DataFrame()
    obj = (mettype==:est) ? lst[1].est : lst[1].obs
    fields = Tuple(Symbol.(prefix .* param_fields(obj)))
    for i ∈ eachindex(lst)
        obj = (mettype==:est) ? lst[i].est : lst[i].obs
        vals = param_values(obj)
        nt = NamedTuple{fields}(vals)
        push!(df, nt)
    end
    df
end


function collect_values(lst::AbstractVector{Met{T}}) where T<:AbstractMetParam
    collect(Iterators.flatten([par.values for par ∈ lst]))
end


function collate_stats(nt; transpose::Bool=false, mettype::Symbol=:est)
    dflst = [collect_mettype(v, "$(k)_", mettype) for (k, v) ∈ pairs(nt)]
    df = hcat(dflst...)
    if transpose
        m = values(nt)
        df.param = string.(getproperty.(m[1], :year))
        df = permutedims(df, "param")
    end
    df
end


function collate_mets(nt)
    mets = values(nt)
    vals = [collect_values(m) for m ∈ mets]

    all_years = map(m -> m.year, mets[1])
    dt = map(enumerate(all_years)) do (i, yr)
        [Date(yr, 1, 1) + Dates.Day(j-1) for j ∈ eachindex(mets[1][i].values)]
    end
    vdt = reduce(vcat, dt)

    parse_date(fn, lst) = (l -> l.value).(fn.(lst))
    year = parse_date(Dates.Year, vdt)
    month = parse_date(Dates.Month, vdt)
    day = parse_date(Dates.Day, vdt)
    doy = dayofyear.(year, month, day)

    nt1 = (; year, month, day, doy)
    nt2 = (; zip(keys(nt), vals)...)
    merged_nt = merge(nt1, nt2)

    DataFrame(merged_nt)
end


function generate_mets(df::AbstractDataFrame; verbose::Bool=true)
    colnames = names(df)

    nt_tmin = (;)
    if any(occursin.("tmin", colnames))
        verbose && println("\nGenerating Tmin")
        tmins = create_temp(df, "tmin")
        foreach(t -> generate!(t; verbose=verbose), tmins)
        nt_tmin = (; tmin=tmins)
    end

    nt_tmax = (;)
    if any(occursin.("tmax", colnames))
        verbose && println("\nGenerating Tmax")
        tmaxs = create_temp(df, "tmax")
        foreach(t -> generate!(t; verbose=verbose), tmaxs)
        nt_tmax = (; tmax=tmaxs)
    end

    nt_wind = (;)
    if any(occursin.("wind", colnames))
        verbose && println("\nGenerating Wind")
        winds = create_wind(df)
        foreach(t -> generate!(t; verbose=verbose), winds)
        nt_wind = (; wind=winds)
    end

    nt_rain = (;)
    if any(occursin.("rain", colnames))
        verbose && println("\nGenerating Rain")
        rains = create_rain(df)
        foreach(t -> generate!(t; verbose=verbose), rains)
        nt_rain = (; rain=rains)
    end

    merge(nt_tmin, nt_tmax, nt_wind, nt_rain)
end


function gof(nt;
             metric::AbstractVector{Function}=[GoF.NMAE, GoF.NMBE, GoF.KGE, GoF.dr],
             ndecimals::Int=-1)
    est = collate_stats(nt; transpose=false, mettype=:est)
    obs = collate_stats(nt; transpose=false, mettype=:obs)

    colnames = propertynames(obs)
    fit = DataFrame(param = string.(colnames))

    vals = Vector{Float64}[]
    for col ∈ colnames
        push!(vals, (g -> g(obs[!, col], est[!, col])).(metric))
    end

    for (i, g) ∈ enumerate(metric)
        fit[!, Symbol(g)] = [v[i] for v ∈ vals]
    end

    if ndecimals > -1
        transform!(fit, Not(1) .=> v -> round.(v; digits=ndecimals), renamecols=false)
    end

    fit
end
