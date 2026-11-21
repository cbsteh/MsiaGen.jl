# The data folder (one folder per site, as MsiaGen reads it), running the generator and
# checking its fit. All the work is done by MsiaGen itself.

const DATA = Ref(normpath(get(ENV, "MSIAGEN_DATA", joinpath(@__DIR__, "..", "..", "data"))))
datadir() = DATA[]
# (a folder name only: never a path that could leave the data folder)
sitedir(s) = (occursin(r"[/\\]", s) || s ∈ (".", "..")) ? throw(ArgumentError("not a site folder: $s")) :
             joinpath(DATA[], s)
obsfile(s) = joinpath(sitedir(s), "$(s)-obs.csv")
statsfile(s) = joinpath(sitedir(s), "$(s)-stats.csv")
simfile(s) = joinpath(sitedir(s), "$(s)-sim.csv")
figfile(s, kind) = joinpath(sitedir(s), "$(s)-$(kind).png")     # kind: "sim" or "fit"

sites() = isdir(DATA[]) ? sort([d for d ∈ readdir(DATA[]) if !startswith(d, ".") && isdir(joinpath(DATA[], d))]) : String[]

# a site's folder name: the name typed, spaces as hyphens (the name MsiaGen's files carry)
site_folder(name) = replace(strip(name), r"\s+" => "-")
# letters, digits, "_", "-" and "."; "." neither first (a hidden folder) nor last (Windows
# drops it), and not a name Windows keeps for a device (CON, NUL, COM1, …, with any extension)
const WIN_RESERVED = r"^(con|prn|aux|nul|com[0-9]|lpt[0-9])(\..*)?$"i
valid_site(name) = (f = site_folder(name);
                    occursin(r"^[A-Za-z0-9_]([A-Za-z0-9_.-]*[A-Za-z0-9_-])?$", f) && !occursin(WIN_RESERVED, f))
const SITE_NAME_RULE = "A site name may hold letters, digits, spaces, '-', '_' and '.' only, may not start or end " *
                       "with '.', and may not be a name Windows keeps for devices (CON, NUL, COM1 and the like)."
# the folder of the site named `name`, or "" if the name is not a valid one (the client sends
# the name, so nothing outside the data folder is read or written)
site_of(name) = valid_site(name) ? site_folder(name) : ""

# latitude line and the years of a file with one (obs, stats), or of the sim file (none)
function file_facts(f; lat_line=true)
    isfile(f) || return nothing
    try
        lat = nothing
        df = open(f) do io
            lat_line && (lat = MsiaGen.parse_lat(readline(io)))
            DataFrame(CSV.File(io; comment="#", ignoreemptyrows=true, select=[:year]))
        end
        ys = sort(unique(df.year))
        (lat=lat, first=first(ys), last=last(ys), n=length(ys))
    catch
        (lat=nothing, first=0, last=0, n=0)
    end
end

span(x) = isnothing(x) ? "" : x.n == 0 ? "unreadable" : x.first == x.last ? "$(x.first)" : "$(x.first)–$(x.last)"

"What a site's folder holds, for its card on the Sites page."
function site_card(s)
    o, st, sm = file_facts(obsfile(s)), file_facts(statsfile(s)), file_facts(simfile(s); lat_line=false)
    lat = something(isnothing(o) ? nothing : o.lat, isnothing(st) ? nothing : st.lat, Some(nothing))
    Dict{String, Any}("site" => s, "lat" => isnothing(lat) ? "" : string(round(lat; digits=3)),
        "obs" => span(o), "stats" => span(st), "sim" => span(sm),
        "fit" => isfile(joinpath(sitedir(s), "$(s)-fit.txt")))
end

"Delete site `s` permanently: its folder and every file in it."
function delete_site(s)
    (valid_site(s) && s in sites()) || error("$s is not a site in the data folder")
    rm(sitedir(s); recursive=true)
end

# ---- generating -----------------------------------------------------------------------------

"""
Run MsiaGen for site `s`: `source` "stats" generates from `<s>-stats.csv` as it is, "obs"
first rebuilds it from `<s>-obs.csv`. Returns the messages MsiaGen printed and the
simulated weather.
"""
function generate(s; source="stats", seed=-1, wetdry_tmax=MsiaGen.WETDRY_TMAX)
    # MsiaGen's messages go to a buffer, not stdout: redirecting stdout (process-wide) from a
    # worker thread would also catch the web server's own output
    io = IOBuffer()
    df = MsiaGen.generate_weather(s; folder=DATA[], seed=seed, from_stats=(source == "stats"),
                                  verbose=false, wetdry_tmax=wetdry_tmax, io=io)
    (log=String(take!(io)), df=df)
end

read_sim(s) = isfile(simfile(s)) ? DataFrame(CSV.File(simfile(s))) : nothing

"One row per year of the simulated weather: annual means and totals."
function sim_summary(df)
    map(collect(groupby(df, :year))) do g
        Dict{String, Any}("year" => g.year[1],
            "tmax" => "tmax" in names(g) ? r2(mean(g.tmax)) : "", "tmin" => "tmin" in names(g) ? r2(mean(g.tmin)) : "",
            "wind" => "wind" in names(g) ? r2(mean(g.wind)) : "",
            "rain" => "rain" in names(g) ? round(Int, sum(g.rain)) : "", "raindays" => "rain" in names(g) ? count(>(0), g.rain) : "")
    end
end

# ---- checking the fit -------------------------------------------------------------------------

"""
How a statistic's fit reads, from the number of its values outside the generator's fit
tolerance (`n_out` of `n`): "good" (under 10% outside), "fair" (under 20%), "acceptable" (30% or less) or "poor" (over 30%).
Shown as a mark in the fit table.
"""
function fit_grade(n_out, n)
    # MsiaGen keeps the best of its attempts for a month that misses, so some misses are normal
    share = n == 0 ? 0.0 : n_out / n
    share < 0.10 ? "good" : share < 0.20 ? "fair" : share <= 0.30 ? "acceptable" : "poor"
end

fmt3(x) = (y = round(x; sigdigits=3); isinteger(y) && abs(y) < 1e6 ? string(Int(y)) : string(y))

"""
The fit of the simulated weather to the stats file of site `s`, as `check_fit` computes it
(without writing its files): the summary rows, one chart per statistic, and what the
simulated weather lacks.
"""
function fit_check(s)
    spec = MsiaGen.csv2df(statsfile(s)).df
    sim = MsiaGen.weather_stats(DataFrame(CSV.File(simfile(s))))
    pairs, gaps = MsiaGen.fit_pairs(spec, sim)
    summ = MsiaGen.fit_summary(pairs)
    pairs.group = [st == "totrain" && m == 0 ? "totrain0" : st for (st, m) ∈ zip(pairs.stat, pairs.month)]
    rows = map(eachrow(summ)) do r
        u = isempty(r.unit) ? "" : " " * r.unit
        Dict{String, Any}("label" => r.label, "unit" => r.unit,
            "mae" => fmt3(r.mae) * u * (ismissing(r.mae_pct) ? "" : " ($(round(r.mae_pct; digits=1))%)"),
            "bias" => (r.bias >= 0 ? "+" : "−") * fmt3(abs(r.bias)) * u,
            "worst" => (r.worst >= 0 ? "+" : "−") * fmt3(abs(r.worst)) * u, "when" => r.worst_when,
            "out" => "$(r.n_out) / $(r.n)", "grade" => fit_grade(r.n_out, r.n))
    end
    charts = map(eachrow(summ)) do r
        g = pairs[pairs.group .== r.group, :]
        st = g.stat[1]
        tolf = v -> MsiaGen.stat_tol(st, v)
        # axis span: at least 10 × the fit tolerance (as MsiaGen's fit chart), but only 4 ×
        # (about 10% of the total) for the annual rain, whose few years would otherwise sit
        # on one spot
        span = r.group == "totrain0" ? 4 * maximum(tolf.(g.spec)) : maximum(MsiaGen.min_span.(st, g.spec))
        fit_chart(g, r.label, r.unit, tolf, span)
    end
    notes = String[]
    isempty(gaps.years) || push!(notes, "The simulated weather lacks year(s) $(join(gaps.years, ", ")) of the stats file; they were not compared.")
    isempty(gaps.stats) || push!(notes, "The simulated weather lacks " * join(MsiaGen.stat_label.(gaps.stats), "; ") * "; not compared.")
    (rows=rows, charts=charts, notes=notes, trend=trend_check(pairs))
end

"Write the report files, as main.jl does: `<s>-fit.txt`, `<s>-fit.png` and `<s>-sim.png`."
function write_reports(s)
    MsiaGen.plot_weather(s; folder=DATA[], show=false, io=devnull)
    MsiaGen.check_fit(s; folder=DATA[], show=false, io=devnull)
end

# ---- is the generated weather up to date? ------------------------------------------------------

when(t) = Dates.format(unix2datetime(t) + (now() - now(UTC)), "d u yyyy, HH:MM")

"""
Where the generated weather of site `s` comes from, for a banner on the Explore and Fit check
tabs: "kind" ("none", "ok" or "warn"), "text" (which file, when it was generated, its years)
and "warn" (why it may not match the statistics or the design now). `dirty`: the design was
changed in the app since it was last started or saved.
"""
function sim_status(s, dirty::Bool)
    f = simfile(s)
    isfile(f) || return Dict{String, Any}("kind" => "none", "warn" => String[],
        "text" => "$s has no generated weather yet. Generate it on the Run MsiaGen tab first.")
    ts = mtime(f)
    facts = file_facts(f; lat_line=false)
    text = "Showing $(basename(f)), generated $(when(ts))" * (isnothing(facts) || facts.n == 0 ? "." : ": $(facts.n) year(s), $(span(facts)).")
    warn = String[]
    st = statsfile(s)
    if !isfile(st)
        push!(warn, "$(basename(st)) is missing, so it is not known which statistics this weather was generated from.")
    elseif mtime(st) > ts + 1
        push!(warn, "$(basename(st)) was changed on $(when(mtime(st))), after this weather was generated, so the weather may not match it. Generate again to bring it up to date.")
    end
    dirty && push!(warn, "Your design has changes that have not been generated yet, so this weather does not include them. To see them, generate again from your design.")
    Dict{String, Any}("kind" => isempty(warn) ? "ok" : "warn", "text" => text, "warn" => warn)
end
