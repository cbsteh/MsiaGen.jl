
function collect_temp(year::Int, data)
    monthdays = days_in_each_month(year)
    mth_data = partition_list(data, monthdays)
    pushfirst!(mth_data, data)
    avg = mean.(mth_data)
    sd = max.(0.01, std.(mth_data))
    rlag = acf1.(mth_data)
    skew = skewness.(mth_data)
    Temp(mean=avg, sd=sd, rlag=rlag, skew=skew)
end


function collect_wind(year::Int, data)
    monthdays = days_in_each_month(year)
    mth_data = partition_list(data, monthdays)
    pushfirst!(mth_data, data)
    avg = mean.(mth_data)
    sd = max.(0.01, std.(mth_data))
    rlag = acf1.(mth_data)
    Wind(mean=avg, sd=sd, rlag=rlag)
end


function collect_rain(year::Int, data)
    totrain, pww, pwd = partition_rain(year, data)
    Rain(totrain=totrain, pww=pww, pwd=pwd)
end


function concat_hdr(namelst, txt::AbstractString)
    lst = [["$(f)$(txt)$(i)" for i ∈ 0:12] for f ∈ namelst]
    join(collect(Iterators.flatten(lst)), ",") * ","
end


function concat_met(met::AbstractMetParam)
    lst = [getproperty(met, f) for f ∈ fieldnames(typeof(met))]
    join(collect(Iterators.flatten(lst)), ",") * ","
end


function create_data_file(wthrfname::AbstractString, inputfname::AbstractString)
    wthr = csv2df(wthrfname)
    df = wthr.df
    fields = ["tmin", "tmax", "wind", "rain"] .∈ (names(df),)

    gdf = groupby(df, :year)

    years = Int[]
    tmins = Temp[]
    tmaxs = Temp[]
    winds = Wind[]
    rains = Rain[]

    for g ∈ gdf
        year = g.year[1]
        push!(years, year)
        fields[1] && push!(tmins, collect_temp(year, g.tmin))
        fields[2] && push!(tmaxs, collect_temp(year, g.tmax))
        fields[3] && push!(winds, collect_wind(year, g.wind))
        fields[4] && push!(rains, collect_rain(year, g.rain))
    end

    hdr = "year,"
    if fields[1] || fields[2]
        namelst = fieldnames(Temp)
        hdr *= fields[1] ? concat_hdr(namelst, "_tmin") : ""
        hdr *= fields[2] ? concat_hdr(namelst, "_tmax") : ""
    end

    hdr *= fields[3] ? concat_hdr(fieldnames(Wind), "_wind") : ""
    hdr *= fields[4] ? concat_hdr(fieldnames(Rain), "") : ""
    hdr = hdr[1:end-1] * "\n"

    vars = "$(wthr.lat)\n"

    fullpath_inputfname = joinpath(dirname(wthrfname), inputfname)
    open(fullpath_inputfname, "w") do fout
        write(fout, vars)
        write(fout, hdr)
        for i ∈ eachindex(years)
            row = "$(years[i])," *
                  (fields[1] ? concat_met(tmins[i]) : "") *
                  (fields[2] ? concat_met(tmaxs[i]) : "") *
                  (fields[3] ? concat_met(winds[i]) : "") *
                  (fields[4] ? concat_met(rains[i]) : "")

            row = row[1:end-1] * "\n"
            write(fout, row)
        end
    end

    fullpath_inputfname
end
