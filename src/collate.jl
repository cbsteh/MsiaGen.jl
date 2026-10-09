
# One table of daily weather from the generated variables in `nt` (all with
# the same years, in year order)
function collate_mets(nt)
    dates = reduce(vcat, [Date(m.year, 1, 1) .+ Day.(0:length(m.values)-1)
                          for m ∈ first(values(nt))])
    vals = (; (k => reduce(vcat, [m.values for m ∈ mets]) for (k, mets) ∈ pairs(nt))...)
    DataFrame(merge((year=Dates.year.(dates), month=Dates.month.(dates),
                     day=Dates.day.(dates), doy=dayofyear.(dates)), vals))
end


# Generate the years in order. Each year starts from the previous year's 31
# December when that year was generated just before it; the first year (or
# a year after a gap) starts from its own start value. `wet`: for Tmax, the
# wet days of each year, made `wetdry` (°C) warmer than dry days.
function generate_years!(mets; verbose::Bool=true, wet=nothing, wetdry=0.0)
    prev, prev_year = nothing, nothing
    for t ∈ sort(mets; by=t -> t.year)
        start = (prev_year == t.year - 1) ? prev : nothing
        kw = isnothing(wet) ? (;) : (wet=wet[t.year], wetdry=wetdry)
        generate!(t; verbose=verbose, prev=start, kw...)
        prev, prev_year = t.values[end], t.year
    end
end


# Repairs any day with tmin >= tmax; returns the number of days repaired
function fix_tmin_tmax!(all_tmins, all_tmaxs)
    n = 0
    for (tn, tx) ∈ zip(all_tmins, all_tmaxs)
        tmins, tmaxs = tn.values, tx.values
        for i ∈ eachindex(tmins)
            if tmaxs[i] < tmins[i]
                tmins[i], tmaxs[i] = tmaxs[i], tmins[i]   # swap positions
                n += 1
            elseif tmaxs[i] == tmins[i]
                tmaxs[i] += 0.1     # slightly increase Tmax; cannot Tmax=Tmin
                n += 1
            end
        end
    end
    n
end


# Variables that can be generated, in order: name, label and how their
# yearly statistics are read from the stats table
const GEN_VARS = ((:tmin, "Tmin", df -> create_temp(df, "tmin")),
                  (:tmax, "Tmax", df -> create_temp(df, "tmax")),
                  (:wind, "Wind", df -> create_wind(df)),
                  (:rain, "Rain", df -> create_rain(df)))


# Generate the variables of the stats table `df`. Rain comes first, so
# that Tmax is `wetdry_tmax` (°C) cooler on wet days than on dry days (0
# for no difference). Returned in the order of GEN_VARS.
function generate_mets(df::AbstractDataFrame; verbose::Bool=true,
                       wetdry_tmax::Real=WETDRY_TMAX)
    check_stats(df)
    colnames = names(df)
    present = [v for v ∈ GEN_VARS if any(occursin.(String(v[1]), colnames))]
    gen = Dict{Symbol,Any}()
    wet = nothing
    for (name, label, create) ∈ sort(present; by=v -> v[1] != :rain)   # rain first
        verbose && println("\nGenerating $(label)")
        mets = create(df)
        kw = (name == :tmax && !isnothing(wet)) ? (wet=wet, wetdry=wetdry_tmax) : (;)
        generate_years!(mets; verbose=verbose, kw...)
        name == :rain && (wet = Dict(m.year => m.values .> 0 for m ∈ mets))
        gen[name] = mets
    end
    nt = NamedTuple{Tuple(first.(present))}(Tuple(gen[first(v)] for v ∈ present))

    # check and repair for any tmin >= tmax occurences:
    if haskey(nt, :tmin) && haskey(nt, :tmax)
        verbose && println("\n\tVerifying Tmin < Tmax")
        n = fix_tmin_tmax!(nt.tmin, nt.tmax)
        n > 0 && println("Tmin was not below Tmax on $(n) day(s); repaired " *
                         "(check_fit measures the repaired weather)")
    end

    nt
end


# Generated months whose best attempt missed the fit tolerance, worst
# first: variable, year, month and score (the largest error relative to
# its tolerance, so above 1)
function misfits(nt)
    df = DataFrame(variable=String[], year=Int[], month=Int[], score=Float64[])
    for (name, label, _) ∈ GEN_VARS
        haskey(nt, name) || continue
        for m ∈ nt[name], (i, s) ∈ enumerate(m.scores)
            s > 1 && push!(df, (label, m.year, i, s))
        end
    end
    sort!(df, :score; rev=true)
end


function print_misfits(io::IO, nt; maxrows::Int=20)
    df = misfits(nt)
    n = sum(length(m.scores) for mets ∈ values(nt) for m ∈ mets; init=0)
    if isempty(df)
        println(io, "All $(n) generated months are within the fit tolerance.")
        return
    end
    println(io, "$(nrow(df)) of $(n) generated months missed the fit tolerance; the " *
                "best attempt was kept (score = largest error / tolerance):")
    for r ∈ first(eachrow(df), maxrows)
        @printf(io, "  %-5s %d %s  %5.2f\n", r.variable, r.year, MONTHS[r.month], r.score)
    end
    nrow(df) > maxrows && println(io, "  … and $(nrow(df) - maxrows) more")
end
