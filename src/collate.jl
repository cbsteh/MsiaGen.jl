
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
# a year after a gap) starts from its own start value.
function generate_years!(mets; verbose::Bool=true)
    prev, prev_year = nothing, nothing
    for t ∈ sort(mets; by=t -> t.year)
        start = (prev_year == t.year - 1) ? prev : nothing
        generate!(t; verbose=verbose, prev=start)
        prev, prev_year = t.values[end], t.year
    end
end


# Repairs any day with tmin >= tmax
function fix_tmin_tmax!(all_tmins, all_tmaxs)
    for (tn, tx) ∈ zip(all_tmins, all_tmaxs)
        tmins, tmaxs = tn.values, tx.values
        for i ∈ eachindex(tmins)
            if tmaxs[i] < tmins[i]
                tmins[i], tmaxs[i] = tmaxs[i], tmins[i]   # swap positions
            elseif tmaxs[i] == tmins[i]
                tmaxs[i] += 0.1     # slightly increase Tmax; cannot Tmax=Tmin
            end
        end
    end
end


# Variables that can be generated, in order: name, label and how their
# yearly statistics are read from the stats table
const GEN_VARS = ((:tmin, "Tmin", df -> create_temp(df, "tmin")),
                  (:tmax, "Tmax", df -> create_temp(df, "tmax")),
                  (:wind, "Wind", df -> create_wind(df)),
                  (:rain, "Rain", df -> create_rain(df)))


function generate_mets(df::AbstractDataFrame; verbose::Bool=true)
    colnames = names(df)
    nt = (;)
    for (name, label, create) ∈ GEN_VARS
        any(occursin.(String(name), colnames)) || continue
        verbose && println("\nGenerating $(label)")
        mets = create(df)
        generate_years!(mets; verbose=verbose)
        nt = merge(nt, NamedTuple{(name,)}((mets,)))
    end

    # check and repair for any tmin >= tmax occurences:
    if haskey(nt, :tmin) && haskey(nt, :tmax)
        verbose && println("\n\tVerifying Tmin < Tmax")
        fix_tmin_tmax!(nt.tmin, nt.tmax)
    end

    nt
end
