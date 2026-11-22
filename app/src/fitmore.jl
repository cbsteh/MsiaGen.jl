# Two more checks on the Fit check tab, beyond the year-by-year comparison:
#
# Trend check: did the generated weather follow the stats file's change over the years?
#   One generated year can miss its value by chance, so for each statistic and month the
#   difference (generated - specified) is fitted with a straight line over the years. Its
#   change from the first to the last year is the trend error: 0 when the generated weather
#   keeps up with the specified change, whatever path that change takes. A month is "off"
#   when the trend error is both beyond chance (outside its 95% range, from the t
#   distribution, so a short record is not judged harshly) and larger than the generator's
#   fit tolerance.
#
# Goal check: counts and lengths that the stats file does not hold, but follow from it
#   (days above a threshold, rainy days, rain per rainy day, dry- and wet-spell lengths):
#   expected from each year's statistics, against the generated daily weather.

# straight-line fit of y on x: slope and its standard error (NaN if too few points)
function line_fit(x, y)
    n = length(x)
    n < 3 && return (slope=NaN, se=NaN)
    xm, ym = mean(x), mean(y)
    sxx = sum((x .- xm) .^ 2)
    iszero(sxx) && return (slope=NaN, se=NaN)
    b = sum((x .- xm) .* (y .- ym)) / sxx
    res = y .- (ym .+ b .* (x .- xm))
    (slope=b, se=sqrt(sum(res .^ 2) / (n - 2) / sxx))
end

const TREND_STATS = ["mean_tmin", "mean_tmax", "sd_tmin", "sd_tmax", "rlag_tmin", "rlag_tmax",
                     "skew_tmin", "skew_tmax", "mean_wind", "sd_wind", "rlag_wind", "totrain", "pww", "pwd"]

signs(x; digits=3) = (y = round(round(x; digits=4); sigdigits=digits); y == 0 ? "0" : (y > 0 ? "+" : "−") * fmt3(abs(y)))

"""
The trend check from the fit pairs (stat, month, year, spec, sim): one row per statistic
and one chart per statistic (designed and generated change of each month). Nothing with
fewer than 5 years.
"""
function trend_check(pairs)
    years = sort(unique(pairs.year))
    length(years) < 5 && return (rows=Dict{String, Any}[], charts=Dict{String, Any}(),
                                 note="The trend check needs at least 5 years: with fewer, a trend cannot be told from chance.")
    span = last(years) - first(years)
    tq = quantile(TDist(length(years) - 2), 0.975)   # 95% range: ± tq standard errors
    rows, charts = Dict{String, Any}[], Dict{String, Any}()
    for st ∈ TREND_STATS
        sel = pairs[(pairs.stat .== st) .& (pairs.month .>= 1), :]
        isempty(sel) && continue
        unit = MsiaGen.stat_unit(st)
        per = map(1:12) do m
            g = sort(sel[sel.month .== m, :], :year)
            x = Float64.(g.year)
            dsg = line_fit(x, g.spec)
            err = line_fit(x, g.sim .- g.spec)
            tol = MsiaGen.stat_tol(st, mean(g.spec))
            e = err.slope * span
            (month=m, designed=dsg.slope * span, generated=(dsg.slope + err.slope) * span, err=e,
             band=tq * err.se * span, off=abs(e) > max(tq * err.se * span, tol))
        end
        any(p -> isnan(p.err), per) && continue
        worst = per[argmax([abs(p.err) for p ∈ per])]
        nout = count(p -> p.off, per)
        u = isempty(unit) ? "" : " " * unit
        push!(rows, Dict{String, Any}("stat" => st, "label" => MsiaGen.stat_label(st),
            "designed" => signs(mean(p.designed for p ∈ per)) * u,
            "generated" => signs(mean(p.generated for p ∈ per)) * u,
            "worst" => signs(worst.err) * u, "when" => MONTH_NAMES[worst.month],
            "out" => "$(nout) / 12", "grade" => fit_grade(nout, 12)))
        charts[st] = trend_chart(per, MsiaGen.stat_label(st), unit, first(years), last(years))
    end
    (rows=rows, charts=charts, note="")
end

# ---- goal check -----------------------------------------------------------------------------

# the stats file's monthly values of one year as named vectors
yrow(df, i, p) = [Float64(df[i, "$(p)$(m)"]) for m ∈ 1:12]

"""
Expected and generated values of each year: days with tmax above `hot` and tmin above
`warm`, rainy days, rain per rainy day, and average dry- and wet-spell lengths.
"""
function goal_values(spec, sim, hot, warm)
    years = sort(intersect(spec.year, sim.year))
    rowof = Dict(y => i for (i, y) ∈ enumerate(spec.year))
    meas = Dict(k => (exp=Float64[], gen=Float64[]) for k ∈ ("hot", "warm", "days", "mm", "dry", "wet"))
    add(k, e, g) = (push!(meas[k].exp, e); push!(meas[k].gen, g))
    for y ∈ years
        i = rowof[y]
        cal = calendar(y)
        d = sim[sim.year .== y, :]
        d = d[sortperm(collect(zip(d.month, d.day))), :]
        for (k, v, thr) ∈ (("hot", "tmax", hot), ("warm", "tmin", warm))
            "mean_$(v)1" ∈ names(spec) && v ∈ names(d) || continue
            μ, σ, γ = yrow(spec, i, "mean_$(v)"), yrow(spec, i, "sd_$(v)"), yrow(spec, i, "skew_$(v)")
            add(k, sum(cal.mdays[m] * p_above(μ[m], σ[m], γ[m], thr) for m ∈ 1:12), count(d[!, v] .> thr))
        end
        ("totrain1" ∈ names(spec) && "rain" ∈ names(d)) || continue
        tot, pww, pwd = yrow(spec, i, "totrain"), yrow(spec, i, "pww"), yrow(spec, i, "pwd")
        # expected: MsiaGen's rain days of each month; spells from the wet/dry chain
        # (dry spells end on a dry-to-wet day, wet spells on a wet-to-dry day)
        W = Float64.(rain_days.(cal.mdays, pww, pwd, tot))
        D = cal.mdays .- W
        wet = d.rain .> 0
        add("days", sum(W), count(wet))
        add("mm", sum(tot) / max(sum(W), 1), sum(d.rain) / max(count(wet), 1))
        add("dry", sum(D) / max(sum(W .* (1 .- pww)), 1e-9), (r = runs(.!wet); isempty(r) ? NaN : mean(r)))
        add("wet", sum(W) / max(sum(D .* pwd), 1e-9), (r = runs(wet); isempty(r) ? NaN : mean(r)))
    end
    years, meas
end

# how close generated must be to expected (on average over the years) to count as close
const GOAL_TOL = Dict("hot" => (abs=3.0, rel=0.10), "warm" => (abs=3.0, rel=0.10), "days" => (abs=5.0, rel=0.05),
                      "mm" => (abs=0.5, rel=0.05), "dry" => (abs=0.2, rel=0.10), "wet" => (abs=0.1, rel=0.10))

"The goal check: one row and one chart per measure."
function goal_check(s, hot, warm)
    spec = MsiaGen.csv2df(statsfile(s)).df
    sim = DataFrame(CSV.File(simfile(s)))
    years, meas = goal_values(spec, sim, hot, warm)
    labels = Dict("hot" => ("Days with tmax above $(hot) °C", "days"), "warm" => ("Days with tmin above $(warm) °C", "days"),
                  "days" => ("Rainy days a year", "days"), "mm" => ("Rain per rainy day", "mm"),
                  "dry" => ("Average dry spell", "days"), "wet" => ("Average wet spell", "days"))
    rows, charts = Dict{String, Any}[], Any[]
    x = Float64.(years)
    span = isempty(years) ? 0 : last(years) - first(years)
    for k ∈ ("hot", "warm", "days", "mm", "dry", "wet")
        e, g = meas[k]
        ok = .!isnan.(g)
        any(ok) || continue
        label, unit = labels[k]
        em, gm = mean(e[ok]), mean(g[ok])
        t = GOAL_TOL[k]
        tol = max(t.abs, t.rel * abs(em))
        near = abs(gm - em) <= tol
        ch(v) = (f = line_fit(x[ok], v[ok]); isnan(f.slope) ? "–" : signs(round(f.slope * span; digits=2); digits=2))
        push!(rows, Dict{String, Any}("label" => label, "unit" => unit,
            "exp" => fmt3(em), "gen" => fmt3(gm), "diff" => signs(gm - em; digits=2),
            "expch" => ch(e), "gench" => ch(g), "grade" => near ? "good" : "poor",
            "mark" => near ? "close" : "off"))
        push!(charts, goal_chart(years, e, g, label, unit, tol))
    end
    (rows=rows, charts=charts)
end

# ---- links between the variables ---------------------------------------------------------------

# the wet-day Tmax shift of the last run of each site in this session (°C)
const RUN_SHIFT = Dict{String, Float64}()

# the measures of MsiaGen's check_links: key, name, unit, and what MsiaGen sets
const LINK_ROWS = [("wet_tmax", "Tmax, wet days minus dry days", "°C"),
                   ("wet_tmin", "Tmin, wet days minus dry days", "°C"),
                   ("r_tmin_tmax", "Tmin–Tmax correlation (r)", ""),
                   ("sd_dtr", "Spread of the daily range, sd(Tmax − Tmin)", "°C")]
const LINK_TOL = 0.15     # °C: how close the generated Tmax wet − dry must be to the shift set

"""
How the variables move together (MsiaGen's `weather_links`, within each month of each year,
averaged over the years): in the generated weather, and in the observed weather when the
site has it. MsiaGen generates each variable on its own, except that Tmax is made `set` °C
warmer on wet days than on dry days (negative: cooler). `known`: whether `set` is the shift
the weather was generated with (or only the value set now). Nothing without tmin, tmax and
rain in the generated weather.
"""
function links_check(s, set, known)
    sim = DataFrame(CSV.File(simfile(s)))
    all(c -> c in names(sim), ("tmin", "tmax", "rain")) || return nothing
    gen = MsiaGen.weather_links(sim)
    obs = nothing
    if isfile(obsfile(s))
        obs = try MsiaGen.weather_links(MsiaGen.csv2df(obsfile(s)).df) catch; nothing end
    end
    f2(x) = isnan(x) ? "–" : replace(string(round(x; digits=2)), "-" => "−")
    rows = map(LINK_ROWS) do (k, label, unit)
        g = gen[end, k]
        o = isnothing(obs) ? NaN : obs[end, k]
        target = k == "wet_tmax" ? f2(set) * " °C" * (known ? "" : " (set now)") :
                 k == "wet_tmin" ? "none (about 0)" : k == "r_tmin_tmax" ? "none (about 0)" : "none"
        near = k == "wet_tmax" ? abs(g - set) <= LINK_TOL : nothing
        Dict{String, Any}("label" => label, "unit" => unit, "obs" => f2(o), "gen" => f2(g), "set" => target,
            "grade" => isnothing(near) ? "" : near ? "good" : "poor", "mark" => isnothing(near) ? "" : near ? "close" : "off")
    end
    # months as 1-12 with one-letter labels, so the 12 fit across a small chart
    x = collect(1:12)
    mtext = MONTH_NAMES
    charts = map(LINK_ROWS) do (k, label, unit)
        tr = Dict{String, Any}[]
        isnothing(obs) || push!(tr, Dict{String, Any}("x" => x, "y" => obs[1:12, k], "mode" => "lines", "name" => "observed",
            "line" => Dict("color" => INK3, "width" => 2), "text" => mtext, "hovertemplate" => "%{text}: observed %{y:.2f}<extra></extra>"))
        k == "wet_tmax" && push!(tr, Dict{String, Any}("x" => x, "y" => fill(set, 12), "mode" => "lines", "name" => "set",
            "line" => Dict("color" => OCHRE, "width" => 1.5, "dash" => "dash"), "hoverinfo" => "skip"))
        push!(tr, Dict{String, Any}("x" => x, "y" => gen[1:12, k], "mode" => "markers", "name" => "generated",
            "marker" => Dict("color" => PRIME, "size" => 7), "text" => mtext, "hovertemplate" => "%{text}: generated %{y:.2f}<extra></extra>"))
        Dict("title" => label, "data" => tr,
             "layout" => base_layout("height" => 260, "showlegend" => false, "margin" => Dict("l" => 56, "r" => 12, "t" => 8, "b" => 36),
                "xaxis" => axis("showgrid" => false, "tickvals" => x, "ticktext" => first.(MONTH_NAMES, 1), "range" => [0.5, 12.5]),
                "yaxis" => axis("title" => Dict("text" => unit), "zeroline" => true, "zerolinecolor" => RULE)))
    end
    (rows=rows, charts=charts, has_obs=!isnothing(obs))
end
