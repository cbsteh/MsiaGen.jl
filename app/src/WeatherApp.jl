module WeatherApp

using Stipple, Stipple.ReactiveTools, StippleUI, StipplePlotly
using CSV, DataFrames, Dates, Random, Sockets, Statistics
using Distributions, SpecialFunctions
import MsiaGen
const Genie = Stipple.Genie
# Vue's <transition> (the step bar's sliding site group) is not an HTML element Genie knows;
# without this it warns UndefVarError(:transition) on the first page load
Genie.Renderer.Html.register_normal_element("transition"; context=Genie.Renderer.Html)

include("designer.jl")
include("charts.jl")
include("sites.jl")
include("simweather.jl")
include("fitmore.jl")

# ---- the controls of each page (label, unit, range, step), drawn by the page template -------

ctl(key, label, unit, min, max, step; hint="", more=false, tkey="") =
    Dict{String, Any}("key" => key, "label" => label, "unit" => unit, "min" => min, "max" => max,
                      "step" => step, "hint" => hint, "more" => more, "tkey" => tkey)

const SPECS = Dict{String, Any}(
    "temp" => Dict(
        "whole" => [ctl("M", "Annual mean", "°C", 15, 40, 0.1; hint="the year's average"),
                    ctl("A", "Amplitude", "°C", 0, 5, 0.05; hint="half the seasonal swing of the monthly means"),
                    ctl("peak", "Warmest day", "day of year", 1, 365, 1; hint="peak of the seasonal curve"),
                    ctl("sd", "sd", "°C", 0.1, 4, 0.05; hint="spread of days around the monthly mean (the band is ±1 sd)"),
                    ctl("rlag", "rlag", "", 0, 0.95, 0.01; hint="how much a day follows the day before: high = smooth, low = jumpy"),
                    ctl("skew", "skew", "", -2, 2, 0.05; hint="negative = more cool dips, positive = more hot spikes")],
        "month" => [ctl("dm", "mean", "°C", -4, 4, 0.1), ctl("dsd", "sd", "°C", -2, 2, 0.05),
                    ctl("drl", "rlag", "", -1, 1, 0.1), ctl("dsk", "skew", "", -2, 2, 0.05)],
        "target" => [ctl("M", "Annual mean", "°C", -8, 8, 0.1), ctl("A", "Amplitude", "°C", -3, 3, 0.05),
                     ctl("peak", "Warmest day", "days", -182, 182, 1),
                     ctl("sd", "sd", "°C", -2, 2, 0.05; more=true), ctl("rlag", "rlag", "", -0.5, 0.5, 0.01; more=true),
                     ctl("skew", "skew", "", -2, 2, 0.05; more=true)],
        # monthly Δ: one row per monthly statistic; tkey is its whole-year change
        "mdelta" => [ctl("dm", "mean", "°C", -8, 8, 0.1; tkey="M"), ctl("dsd", "sd", "°C", -2, 2, 0.05; tkey="sd"),
                     ctl("drl", "rlag", "", -0.5, 0.5, 0.01; tkey="rlag"), ctl("dsk", "skew", "", -2, 2, 0.05; tkey="skew")],
        "rcols" => ["mean (°C)", "sd (°C)", "rlag", "skew"]),
    "wind" => Dict(
        "whole" => [ctl("M", "Annual mean", "m/s", 0.1, 10, 0.05; hint="the year's average"),
                    ctl("A", "Amplitude", "m/s", 0, 3, 0.05; hint="half the seasonal swing of the monthly means"),
                    ctl("peak", "Windiest day", "day of year", 1, 365, 1; hint="peak of the seasonal curve; the northeast monsoon is the windiest time"),
                    ctl("sd", "sd", "m/s", 0.05, 2, 0.01; hint="spread of days; sd / mean sets how gusty the spikes are"),
                    ctl("rlag", "rlag", "", 0, 0.95, 0.01; hint="how much a day follows the day before")],
        "month" => [ctl("dm", "mean", "m/s", -2, 2, 0.05), ctl("dsd", "sd", "m/s", -1, 1, 0.01),
                    ctl("drl", "rlag", "", -1, 1, 0.1)],
        "target" => [ctl("M", "Annual mean", "m/s", -3, 3, 0.05), ctl("A", "Amplitude", "m/s", -1.5, 1.5, 0.05),
                     ctl("peak", "Windiest day", "days", -182, 182, 1),
                     ctl("sd", "sd", "m/s", -1, 1, 0.01; more=true), ctl("rlag", "rlag", "", -0.5, 0.5, 0.01; more=true)],
        "mdelta" => [ctl("dm", "mean", "m/s", -3, 3, 0.05; tkey="M"), ctl("dsd", "sd", "m/s", -1, 1, 0.01; tkey="sd"),
                     ctl("drl", "rlag", "", -0.5, 0.5, 0.01; tkey="rlag")],
        "rcols" => ["mean (m/s)", "sd (m/s)", "rlag"]),
    "rain" => Dict(
        "whole" => [ctl("total", "Annual rainfall", "mm", 100, 6000, 10; hint="the year's total"),
                    ctl("seas", "Seasonality", "%", 0, 90, 5; hint="at 30%, the driest months get about 0.7× and the wettest 1.3× the average month"),
                    ctl("dry", "Driest day", "day of year", 1, 365, 1; hint="middle of the dry season: mid-year on the west coast, day 45–75 in the south and east"),
                    ctl("pww", "pww", "", 0.01, 0.99, 0.01; hint="chance of a wet day after a wet day: higher = longer wet spells"),
                    ctl("pwd", "pwd", "", 0.01, 0.99, 0.01; hint="chance of a wet day after a dry day: lower = longer dry spells")],
        "month" => [ctl("dr", "rainfall", "%", -100, 200, 1), ctl("dww", "pww", "", -0.5, 0.5, 0.01),
                    ctl("dwd", "pwd", "", -0.5, 0.5, 0.01)],
        "target" => [ctl("total", "Annual rainfall", "%", -80, 100, 1), ctl("seas", "Seasonality", "% points", -50, 50, 5),
                     ctl("dry", "Driest day", "days", -182, 182, 1),
                     ctl("pww", "pww", "", -0.5, 0.5, 0.01; more=true), ctl("pwd", "pwd", "", -0.5, 0.5, 0.01; more=true)],
        "mdelta" => [ctl("dr", "rainfall", "%", -80, 100, 1; tkey="total"), ctl("dww", "pww", "", -0.5, 0.5, 0.01; tkey="pww"),
                     ctl("dwd", "pwd", "", -0.5, 0.5, 0.01; tkey="pwd")],
        "rcols" => ["rainfall (mm)", "pww", "pwd", "rain days", "wet spell (d)", "dry spell (d)"]))

const PATH_OPTS = [Dict("label" => l, "value" => v) for (v, l) ∈ PATHS]

# ---- the design and its baseline ------------------------------------------------------------

"The design of all four parts, from a file's average year `b` (or typical values: nothing)."
function new_design(b)
    get_b(k) = b isa NamedTuple ? getproperty(b, k) : nothing
    Dict{String, Any}("tmin" => temp_design("tmin", get_b(:tmin)), "tmax" => temp_design("tmax", get_b(:tmax)),
                      "wind" => wind_design(get_b(:wind)), "rain" => rain_design(get_b(:rain)))
end

const MAX_UPLOAD = 20_000_000       # bytes: a century of daily weather is about 2 MB
const UPLOAD = Ref((name="", text=""))   # the last daily weather file uploaded on the Site page

"The daily weather file a design starts from: (name, text), or nothing for typical values."
function baseline_text(src, site)
    src == "obs" && isfile(obsfile(site)) && return (basename(obsfile(site)), read(obsfile(site), String))
    src == "sim" && isfile(simfile(site)) && return (basename(simfile(site)), read(simfile(site), String))
    src == "upload" && !isempty(UPLOAD[].text) && return (UPLOAD[].name, UPLOAD[].text)
    nothing
end

function base_options(site)
    o = [Dict("label" => "Typical Malaysian values", "value" => "typical")]
    isfile(obsfile(site)) && push!(o, Dict("label" => "Observed weather ($(site)-obs.csv)", "value" => "obs"))
    isfile(simfile(site)) && push!(o, Dict("label" => "Simulated weather ($(site)-sim.csv)", "value" => "sim"))
    isempty(UPLOAD[].text) || push!(o, Dict("label" => "Uploaded file ($(UPLOAD[].name))", "value" => "upload"))
    o
end

"""
Start the design again from the chosen baseline: every slider from the file's average
year (or typical values), and the latitude and years from the file.
"""
function apply_baseline!(m)
    f = baseline_text(m.base_src[], site_of(m.site[]))
    b = nothing
    if isnothing(f)
        m.base_note[] = Dict("kind" => "info", "text" => "No file: every control starts from typical Malaysian values.")
    else
        b = try
            average_year(parse_weather(f[2]))
        catch e
            m.base_note[] = Dict("kind" => "negative", "text" => "Could not read $(f[1]): $(sprint(showerror, e))")
            nothing
        end
        if !isnothing(b)
            vars = ("tmin", "tmax", "wind", "rain")
            found = [v for v ∈ vars if !isnothing(getproperty(b, Symbol(v)))]
            absent = setdiff(collect(vars), found)
            m.base_note[] = Dict("kind" => "positive", "text" =>
                "Average year of $(f[1]): $(length(b.years)) year(s), $(first(b.years))–$(last(b.years)). " *
                "The controls start from it for $(join(found, ", "))" *
                (isempty(absent) ? "." : " ($(join(absent, ", ")) missing or incomplete in the file: typical values used).") *
                (isnothing(b.lat) ? " The file has no latitude line, so check the latitude." : ""))
            isnothing(b.lat) || (m.lat[] = clamp(round(b.lat; digits=3), -10.0, 10.0))
            m.start_year[] = clamp(first(b.years), 1950, 2100)
            m.end_year[] = clamp(last(b.years), 1950, 2100)
        end
    end
    m.hot_goal[] = -1.0
    m.rs_dry[] = -1.0
    m.rs_mm[] = -1.0
    d = new_design(b)
    # separate copies, so an edit of the design cannot also change what it is compared with
    m.design0[] = deepcopy(d)
    DESIGN0[] = design_state(d, m)
    m.design[] = d
    m.design_dirty[] = false
end

# What the stats file is written from: the design, the latitude and the years. DESIGN0 holds
# it as last started or saved, to tell if the design was changed since.
design_state(d, m) = (design=deepcopy(d), lat=m.lat[], years=(m.start_year[], m.end_year[]))
const DESIGN0 = Ref(design_state(new_design(nothing), (lat=Ref(0.0), start_year=Ref(0), end_year=Ref(0))))
is_dirty(m) = (s = DESIGN0[]; m.design[] != s.design || m.lat[] != s.lat ||
                              (m.start_year[], m.end_year[]) != s.years)

# ---- what each page shows -------------------------------------------------------------------

design_years(m) = m.start_year[]:max(m.start_year[], m.end_year[])
signed(v; digits=3) = (v = round(v; digits); v == 0 ? "" : (v > 0 ? "+" : "−") * string(abs(v)))
cell(v, d=0.0; digits=3) = Dict("v" => string(v), "d" => signed(d; digits))
capped_note(ys, part) = isempty(ys) ? String[] :
    ["Some $(part) values went past their limits in $(length(ys)) year(s), $(first(ys))–$(last(ys)). They were capped at the limit."]

# true if the design part `d` changes by the end year (any end-year change not 0)
has_change(d) = any(k != "path" && num(v) != 0 for (k, v) ∈ d["target"]) ||
                any(any(!=(0), num.(v)) for v ∈ values(get(d, "mdelta", Dict())))

"""
The chart and table of a variable page: the start and end year on one chart (each shown
when its box is ticked; the end year only when the design changes by then), and the
monthly values of one of the two years.
"""
function var_page!(m)
    nav = m.nav[]
    nav ∈ ("temp", "wind", "rain") || return
    years = design_years(m)
    st = all_stats(m.design[], years)
    part = nav == "temp" ? m.tvar[] : nav
    changed = length(years) > 1 && has_change(m.design[][part])
    m.vchanged[] = changed
    d = m.design[][part]["month"]
    show_s = m.show_start[]
    show_e = changed && m.show_end[]
    grids = Any[]
    if nav == "temp"
        v = part
        S, E = first(st.t[v]), last(st.t[v])
        dr = t_draws(m.t_new[])
        ps, pe = t_sample_year(S, dr), t_sample_year(E, dr)
        ov = other_tvar(v)
        sel, oth = Any[], Any[]
        if show_s
            push!(sel, (S, ps, YEAR_COLORS[v][1], "$(S.year) start year"))
            m.show_other[] && push!(oth, (first(st.t[ov]), YEAR_COLORS[ov][1], "$(uppercasefirst(ov)) $(S.year)"))
        end
        if show_e
            push!(sel, (E, pe, YEAR_COLORS[v][2], "$(E.year) end year"))
            m.show_other[] && push!(oth, (last(st.t[ov]), YEAR_COLORS[ov][2], "$(uppercasefirst(ov)) $(E.year)"))
        end
        chart = isempty(sel) ? nothing : overlay_chart(sel; unit="°C", ylabel="$(uppercasefirst(v)) (°C)", other=oth, days=m.show_days[])
        rows = (s, ch) -> [Dict("month" => MONTH_NAMES[k], "cells" => [cell(r2(s.μ[k]), ch ? d["dm"][k] : 0.0),
            cell(s.σ[k], ch ? d["dsd"][k] : 0.0), cell(s.ρ[k], ch ? d["drl"][k] : 0.0), cell(s.γ[k], ch ? d["dsk"][k] : 0.0)]) for k ∈ 1:12]
        summ = s -> "Whole year $(s.year): annual mean $(r2(s.M)) °C, amplitude $(r2(s.A)) °C, warmest day $(show_day(s.peak))."
        thr = m.hot_thr[]
        m.hot_days[] = Dict("start" => round(days_above(S, thr); digits=1), "end" => round(days_above(E, thr); digits=1))
        n = t_cross_note(st.t)
        m.vnotes[] = [capped_note(capped_years(st).temperature, "temperature"); isnothing(n) ? String[] : [n]]
    elseif nav == "wind"
        S, E = first(st.w), last(st.w)
        dr = w_draws(m.w_new[])
        ps, pe = w_sample_year(S, dr), w_sample_year(E, dr)
        sel = Any[]
        show_s && push!(sel, (S, ps, C_START, "$(S.year) start year"))
        show_e && push!(sel, (E, pe, C_END, "$(E.year) end year"))
        chart = isempty(sel) ? nothing : overlay_chart(sel; unit="m/s", ylabel="Wind speed (m/s)", floor_at=W_FLOOR,
                                                       floor_label="0.1 m/s floor", days=m.show_days[])
        rows = (s, ch) -> [Dict("month" => MONTH_NAMES[k], "cells" => [cell(r2(s.μ[k]), ch ? d["dm"][k] : 0.0),
            cell(s.σ[k], ch ? d["dsd"][k] : 0.0), cell(s.ρ[k], ch ? d["drl"][k] : 0.0)]) for k ∈ 1:12]
        summ = s -> "Whole year $(s.year): annual mean $(r2(s.M)) m/s, amplitude $(r2(s.A)) m/s, windiest day $(show_day(s.peak))."
        n = w_floor_note(st.w)
        m.vnotes[] = [capped_note(capped_years(st).wind, "wind"); isnothing(n) ? String[] : [n]]
    else
        S, E = first(st.r), last(st.r)
        dr = r_draws(m.r_new[])
        sel = Any[]
        show_s && push!(sel, (S, C_START, "$(S.year) start year", 1))
        show_e && push!(sel, (E, C_END, "$(E.year) end year", 2))
        chart = isempty(sel) ? nothing : rain_overlay_chart(sel)
        grids = [Dict{String, Any}("label" => l, "color" => c, "chart" => rain_grid_chart(y, make_strip(y, dr))) for (y, c, l, _) ∈ sel]
        rows = (s, ch) -> [Dict("month" => MONTH_NAMES[k], "cells" => [cell(round(Int, s.total[k]), ch ? d["dr"][k] : 0.0; digits=0),
            cell(s.pww[k], ch ? d["dww"][k] : 0.0), cell(s.pwd[k], ch ? d["dwd"][k] : 0.0), cell(s.raindays[k]),
            cell(round(s.wetspell[k]; digits=1)), cell(round(s.dryspell[k]; digits=1))]) for k ∈ 1:12]
        summ = s -> "Whole year $(s.year): annual rainfall $(round(Int, sum(s.total))) mm, " *
                  "$(sum(s.raindays)) rain days, seasonality $(round(Int, s.seas))%, driest day $(show_day(s.dry))."
        m.vnotes[] = capped_note(capped_years(st).rain, "rainfall")
    end
    # the table: the start year with each month's change, or the end year
    show_end = changed && m.view_end[]
    m.vrows[] = show_end ? rows(E, false) : rows(S, true)
    m.vsummary[] = summ(show_end ? E : S)
    m.vchart[] = something(chart, Dict("data" => [], "layout" => Dict()))
    m.vgrids[] = grids
end

stats_name(site) = "$(site_folder(site))-stats.csv"

function stats_page!(m)
    m.nav[] == "stats" || return
    years = design_years(m)
    st = all_stats(m.design[], years)
    m.checks[] = value_checks(st)
    annual(s) = r2(sum(s.μ .* s.mdays) / s.ndays)
    m.annual_rows[] = [Dict("year" => y, "tmax" => annual(st.t["tmax"][i]), "tmin" => annual(st.t["tmin"][i]),
                            "wind" => annual(st.w[i]), "rain" => round(Int, sum(st.r[i].total)), "raindays" => sum(st.r[i].raindays))
                       for (i, y) ∈ enumerate(years)]
    m.stats_text[] = stats_csv(m.lat[], years, st)
    m.trends[] = trend_charts(st, years)
    s = site_of(m.site[])
    m.stats_exists[] = !isempty(s) && isfile(statsfile(s))
end

"Problems that stop the stats file being written (empty if none)."
function design_problems(m)
    p = String[]
    isempty(strip(m.site[])) && push!(p, "Give the site a name (Site page).")
    !isempty(strip(m.site[])) && !valid_site(m.site[]) && push!(p, SITE_NAME_RULE)
    m.start_year[] > m.end_year[] && push!(p, "The start year is after the end year (Site page).")
    p
end

"Which files the Generate-from options need are there; an option without its file is not chosen."
function refresh_run_has!(m)
    s = site_of(m.site[])
    h = Dict("stats" => !isempty(s) && isfile(statsfile(s)), "obs" => !isempty(s) && isfile(obsfile(s)))
    m.run_has[] = h
    get(h, m.run_src[], true) || (m.run_src[] = "design")
end

function save_stats!(m)
    s = site_folder(m.site[])
    mkpath(sitedir(s))
    years = design_years(m)
    write(statsfile(s), stats_csv(m.lat[], years, all_stats(m.design[], years)))
    # the saved design is now the stats file, so it no longer counts as changed (otherwise the
    # Explore and Fit check tabs warn that the design was not generated, right after generating it)
    DESIGN0[] = design_state(m.design[], m)
    m.design_dirty[] = false
    m.site_cards[] = site_card.(sites())
    m.stats_exists[] = true
    refresh_run_has!(m)
    s
end

function load_sim!(m, s)
    df = read_sim(s)
    if isnothing(df)
        m.sim_rows[] = Dict{String, Any}[]
        return
    end
    m.sim_rows[] = sim_summary(df)
end

const FIT_TIME = Ref(0.0)   # when the fit shown was computed
mtime_or0(f) = isfile(f) ? mtime(f) : 0.0

function load_fit!(m, s)
    m.fit_site[] = s
    FIT_TIME[] = time()
    if !(isfile(statsfile(s)) && isfile(simfile(s)))
        m.fit_rows[] = Dict{String, Any}[]; m.fit_charts[] = Any[]
        m.fit_notes[] = String[]
        m.trend_rows[] = Dict{String, Any}[]; m.trend_charts[] = Dict{String, Any}(); m.trend_note[] = ""
        m.goal_rows[] = Dict{String, Any}[]; m.goal_charts[] = Any[]
        return
    end
    try
        f = fit_check(s)
        m.fit_rows[] = f.rows; m.fit_charts[] = f.charts; m.fit_notes[] = f.notes
        m.trend_rows[] = f.trend.rows; m.trend_charts[] = f.trend.charts; m.trend_note[] = f.trend.note
    catch e
        m.fit_rows[] = Dict{String, Any}[]; m.fit_charts[] = Any[]
        m.trend_rows[] = Dict{String, Any}[]; m.trend_charts[] = Dict{String, Any}(); m.trend_note[] = ""
        m.fit_notes[] = ["Could not check the fit: $(sprint(showerror, e))"]
    end
    load_goals!(m, s)
    m.figs[] = Dict("sim" => isfile(figfile(s, "sim")), "fit" => isfile(figfile(s, "fit")))
    m.fitrev[] += 1
end

"The goal check of site `s` (empty without a stats file and generated weather)."
function load_goals!(m, s)
    m.goal_rows[] = Dict{String, Any}[]; m.goal_charts[] = Any[]
    isfile(statsfile(s)) && isfile(simfile(s)) || return
    try
        g = goal_check(s, Float64(m.goal_hot[]), Float64(m.goal_warm[]))
        m.goal_rows[] = g.rows; m.goal_charts[] = g.charts
    catch e
        m.fit_notes[] = [m.fit_notes[]; "Could not run the goal check: $(sprint(showerror, e))"]
    end
    m.goalrev[] += 1
end

const SW_REV = Ref(0)

"The Explore tab (the generated weather) of the site open (when that page is shown)."
function sim_page!(m)
    s = site_of(m.site[])
    m.nav[] ∈ ("simw", "fit") && !isempty(s) && (m.sim_info[] = sim_status(s, m.design_dirty[]))
    m.nav[] == "simw" || return
    r =isempty(s) ? nothing : sim_page(s, m.sw_var[], m.sw_year[],
                                        (hot=Float64(m.th_hot[]), dry_day=Float64(m.th_dry_day[]), dry_month=Float64(m.th_dry_month[])))
    if isnothing(r)
        m.sw_years[] = Int[]; m.sw_charts[] = Any[]; m.sw_stats[] = Dict{String, Any}[]; m.sw_ext[] = Dict{String, Any}[]
        m.sw_chances[] = Dict{String, Any}[]
        m.sw_note[] = "$(isempty(s) ? "This site" : s) has no simulated weather yet. Generate it first."
        return
    end
    m.sw_year[] in [0; r.years] || (m.sw_year[] = 0)
    # a fresh key for every chart on each update, so the browser draws new charts rather than
    # updating the old ones in place (Plotly.react can leave a chart's box collapsed, and the
    # charts below then overlap it)
    rev = SW_REV[] += 1
    foreach(c -> c["key"] = "$(rev)-$(c["key"])", r.charts)
    m.sw_years[] = r.years; m.sw_charts[] = r.charts; m.sw_stats[] = r.stats; m.sw_ext[] = r.extremes; m.sw_chances[] = r.chances; m.sw_note[] = r.note
end

"Open site `s`: its name, and a design started from its observed weather if it has any."
function open_site!(m, s)
    m.site[] = s
    m.is_new[] = false
    m.base_opts[] = base_options(s)
    facts = something(file_facts(obsfile(s)), file_facts(statsfile(s)), Some(nothing))
    isnothing(facts) || isnothing(facts.lat) || (m.lat[] = clamp(round(facts.lat; digits=3), -10.0, 10.0))
    isnothing(facts) || facts.n == 0 || (m.start_year[] = facts.first; m.end_year[] = facts.last)
    src = isfile(obsfile(s)) ? "obs" : "typical"
    m.base_src[] == src ? apply_baseline!(m) : (m.base_src[] = src)   # a change of base_src applies it
    m.run_src[] = isfile(statsfile(s)) ? "stats" : "design"
    refresh_run_has!(m)
    load_sim!(m, s)
end

toastmsg(msg, kind="positive") = Dict("msg" => msg, "type" => kind, "t" => time())

const ABOUT = Dict(
    "authors" => "Teh, C. B. S., Cheah, S. S., & Appleton, D. R.", "year" => "2026",
    "title" => "A stochastic daily weather generator for perennial crop simulations in tropical Malaysia",
    "journal" => "PLoS One", "volume" => "21(2)", "article" => "e0338833",
    "doi" => "https://doi.org/10.1371/journal.pone.0338833", "github" => "https://github.com/cbsteh/MsiaGen.jl")

# ---- reactive model --------------------------------------------------------------------------

@app begin
    @in nav = "home"
    @out about = ABOUT
    @out toast = Dict("msg" => "", "type" => "positive", "t" => 0.0)
    @out specs = SPECS
    @out path_opts = PATH_OPTS
    @out data_folder = datadir()
    @out site_cards = site_card.(sites())
    @in open_site = ""
    @out opening = ""
    @in new_dialog = false
    @in new_name = ""
    @in do_new = false
    # deleting a site: the dialog is open while del_target holds its name; step 1 warns and
    # wants the name typed exactly, step 2 asks once more
    @in del_target = ""
    @in del_step = 1
    @in del_typed = ""
    @in do_delete = false
    @out del_files = String[]

    # site and baseline
    @in site = ""
    @out is_new = true
    @in lat = 3.0
    @in start_year = 2000
    @in end_year = 2020
    @in base_src = "typical"
    @out base_opts = [Dict("label" => "Typical Malaysian values", "value" => "typical")]
    @out base_note = Dict("kind" => "info", "text" => "No file: every control starts from typical Malaysian values.")
    @in upload_rev = 0
    @in restart_design = false

    # the design (one Dict per part) and what the variable pages show
    @in design = new_design(nothing)
    @out design_dirty = false
    @out design0 = new_design(nothing)   # the design as it started, for undoing a month's changes
    @in tvar = "tmin"
    @in ctab = "whole"
    @in view_end = false
    @in t_new = 0
    @in w_new = 0
    @in r_new = 0
    @out vchart = Dict("data" => [], "layout" => Dict())
    @out vgrids = Any[]
    @out year_colors = Dict(k => collect(v) for (k, v) ∈ YEAR_COLORS)
    @out rain_legend = [Dict("color" => c, "label" => l) for (c, l) ∈ zip(RAIN_COLORS, RAIN_LABELS)]
    @in show_start = true
    @in show_end = true
    @in show_days = true
    @out vchanged = false
    @in show_other = false
    @out vrows = Dict{String, Any}[]
    @out vsummary = ""
    @out vnotes = String[]
    # goals on the End year tab: days above a threshold (temperature) and the rain pattern;
    # a goal of -1 means "not set" (the box shows the end year as designed)
    @in hot_thr = 25.0          # tmin is shown first; tmax switches it to 35
    @in hot_goal = -1.0
    @in hot_by = "M"
    @in apply_hot = false
    @out hot_days = Dict("start" => 0.0, "end" => 0.0)
    @in rs_dry = -1.0
    @in rs_mm = -1.0
    @in md_key = ""

    # statistics
    @out checks = String[]
    @out annual_rows = Dict{String, Any}[]
    @out stats_text = ""
    @out trends = Any[]
    @out stats_exists = false
    @in save_stats = false
    @in stats_view = "table"

    # generate
    @in run_src = "design"
    @out run_has = Dict("stats" => false, "obs" => false)
    @in seed = -1
    @in wetdry = MsiaGen.WETDRY_TMAX
    @in start_run = false
    @out running = false
    @out run_log = ""
    @out sim_rows = Dict{String, Any}[]

    # simulated weather
    @in sw_var = "tmin"
    @in sw_year = 0
    @out sw_years = Int[]
    @out sw_charts = Any[]
    @out sw_stats = Dict{String, Any}[]
    @out sw_ext = Dict{String, Any}[]
    @out sw_chances = Dict{String, Any}[]
    @out sw_note = ""
    @out sim_info = Dict{String, Any}("kind" => "none", "text" => "", "warn" => String[])
    @in th_hot = 33.0
    @in th_dry_day = 0.5
    @in th_dry_month = 100.0

    # fit check
    @out fit_site = ""
    @out fit_rows = Dict{String, Any}[]
    @out fit_charts = Any[]
    @out fit_notes = String[]
    @out fitrev = 0
    # trend check (one chart per statistic, the one picked shown) and goal check
    @out trend_rows = Dict{String, Any}[]
    @out trend_charts = Dict{String, Any}()
    @out trend_note = ""
    @in trend_stat = "mean_tmax"
    @out goal_rows = Dict{String, Any}[]
    @out goal_charts = Any[]
    @out goalrev = 0
    @in goal_hot = 35.0
    @in goal_warm = 25.0
    @in make_reports = false
    @out reporting = false
    @out figs = Dict("sim" => false, "fit" => false)
    @out figrev = 0

    @onchange design begin
        design_dirty = is_dirty(__model__)
        var_page!(__model__)
        stats_page!(__model__)
    end
    @onchange tvar begin
        hot_goal = -1.0
        thr = tvar == "tmax" ? 35.0 : 25.0
        # a new hot_thr redraws the page through its own handler; redraw here only if it stays
        thr == hot_thr ? var_page!(__model__) : (hot_thr = thr)
    end
    @onchange hot_thr begin
        var_page!(__model__)
    end
    @onbutton apply_hot begin
        years = design_years(__model__)
        if length(years) < 2
            toast = toastmsg("Set an end year after the start year first (on the Site page).", "warning")
        elseif hot_goal >= 0
            c = only(filter(c -> c["key"] == hot_by, SPECS["temp"]["target"]))
            x, n = hot_change(design[tvar], last(years), hot_thr, hot_goal, hot_by, c["min"], c["max"], c["step"])
            nd = deepcopy(design)
            nd[tvar]["target"][hot_by] = x
            what = hot_by == "M" ? "annual mean" : "sd"
            # the goal is out of reach only when the change sits at its slider's end; otherwise
            # the small gap comes from the slider step
            far = x in (c["min"], c["max"]) && abs(n - hot_goal) > 0.5 ?
                  " (the most the $(what) change allows: the goal is out of reach)" :
                  " (the nearest the $(c["step"]) °C steps give)"
            hot_goal = -1.0
            design = nd
            toast = toastmsg("$(uppercasefirst(tvar)): $(what) change set to $(signed(x)) °C, giving about " *
                             "$(round(n; digits=1)) days above $(hot_thr) °C in $(last(years))$(far).")
        end
    end
    @onchange start_year, end_year begin
        design_dirty = is_dirty(__model__)
    end
    @onchange nav, view_end, show_other, show_start, show_end, show_days, t_new, w_new, r_new, start_year, end_year begin
        var_page!(__model__)
        stats_page!(__model__)
        nav == "run" && refresh_run_has!(__model__)
        s = site_of(site)
        if nav == "fit" && !isempty(s)
            fit_site == s && !isempty(fit_rows) && FIT_TIME[] >= mtime_or0(simfile(s)) && FIT_TIME[] >= mtime_or0(statsfile(s)) ||
                load_fit!(__model__, s)
        end
    end
    @onchange goal_hot, goal_warm begin
        nav == "fit" && !isempty(fit_site) && load_goals!(__model__, fit_site)
    end
    @onchange lat, site begin
        design_dirty = is_dirty(__model__)
        stats_page!(__model__)
    end
    @onchange base_src begin
        apply_baseline!(__model__)
        toast = toastmsg("The design starts from " *
                         (base_src == "typical" ? "typical Malaysian values." : "the chosen file's average year."), "info")
    end
    @onbutton restart_design begin
        apply_baseline!(__model__)
        toast = toastmsg("The design was started again.", "info")
    end
    @onchange upload_rev begin
        if upload_rev > 0
            base_opts = base_options(site_of(site))
            base_src == "upload" ? apply_baseline!(__model__) : (base_src = "upload")
        end
    end
    @onchange open_site begin
        if !isempty(open_site)
            ok = open_site in sites()      # only a site of the data folder (the name comes from the browser)
            ok && open_site!(__model__, open_site)
            open_site = ""
            ok && (nav = "site")
        end
    end
    @onbutton do_new begin
        name = site_folder(new_name)
        if !valid_site(name)
            toast = toastmsg(SITE_NAME_RULE, "negative")
        elseif name in sites()
            toast = toastmsg("$name already exists: open it from its card.", "warning")
        else
            site = name
            is_new = true
            base_opts = base_options(name)
            base_src == "typical" ? apply_baseline!(__model__) : (base_src = "typical")
            run_src = "design"
            refresh_run_has!(__model__)
            load_sim!(__model__, name)
            new_dialog = false
            new_name = ""
            nav = "site"
        end
    end
    @onchange del_target begin
        del_files = isempty(del_target) || !(del_target in sites()) ? String[] : sort(readdir(sitedir(del_target)))
    end
    @onbutton do_delete begin
        s = del_target
        if isempty(s) || del_typed != s || del_step != 2
            toast = toastmsg("Type the site name exactly to delete it.", "warning")
        else
            try
                delete_site(s)
                if s == site_folder(site)     # the site open now: back to no site
                    site = ""
                    is_new = true
                    fit_site = ""
                    base_opts = base_options("")
                    base_src == "typical" ? apply_baseline!(__model__) : (base_src = "typical")
                end
                site_cards = site_card.(sites())
                toast = toastmsg("$s and its files were deleted.", "info")
            catch e
                toast = toastmsg("Could not delete $s: $(sprint(showerror, e))", "negative")
            end
            del_target = ""
            del_typed = ""
            del_step = 1
        end
    end
    @onbutton save_stats begin
        p = design_problems(__model__)
        if isempty(p)
            s = save_stats!(__model__)
            is_new = false
            toast = toastmsg("Saved as $(stats_name(s)) in $(sitedir(s)).")
        else
            toast = toastmsg(join(p, " "), "negative")
        end
    end

    @onbutton start_run begin
        s = site_of(site)
        p = run_src == "design" ? design_problems(__model__) :
            isempty(strip(site)) ? ["Open or create a site first."] : isempty(s) ? [SITE_NAME_RULE] : String[]
        run_src == "stats" && !isfile(statsfile(s)) && push!(p, "$s has no $(s)-stats.csv yet.")
        run_src == "obs" && !isfile(obsfile(s)) && push!(p, "$s has no observed weather ($(s)-obs.csv).")
        if !isempty(p)
            toast = toastmsg(join(p, " "), "negative")
        else
            running = true
            run_log = ""
            try
                if run_src == "design"
                    save_stats!(__model__)
                    is_new = false
                end
                sd = Int(seed)
                w = Float64(wetdry)
                src = run_src == "obs" ? "obs" : "stats"
                r = fetch(Threads.@spawn generate(s; source=src, seed=sd, wetdry_tmax=w))
                run_log = r.log
                site_cards = site_card.(sites())
                base_opts = base_options(s)
                load_sim!(__model__, s)
                load_fit!(__model__, s)
                sim_page!(__model__)
                toast = toastmsg("Weather generated. Click Fit check to see how well it matches your design, or Explore to see what it looks like.")
            catch e
                @error "generate" exception = (e, catch_backtrace())
                run_log = sprint(showerror, e)
                toast = toastmsg("MsiaGen stopped: see the messages below.", "negative")
            finally
                running = false
            end
        end
    end
    @onchange nav, site, design_dirty, sw_var, sw_year, th_hot, th_dry_day, th_dry_month begin
        sim_page!(__model__)
    end
    @onbutton make_reports begin
        s = fit_site
        if isfile(statsfile(s)) && isfile(simfile(s))
            reporting = true
            try
                fetch(Threads.@spawn write_reports(s))
                figs = Dict("sim" => isfile(figfile(s, "sim")), "fit" => isfile(figfile(s, "fit")))
                figrev = figrev + 1
                site_cards = site_card.(sites())
                toast = toastmsg("Saved $(s)-fit.txt, $(s)-fit.png and $(s)-sim.png in $(sitedir(s)).")
            catch e
                @error "reports" exception = (e, catch_backtrace())
                toast = toastmsg("Could not write the report files: $(sprint(showerror, e))", "negative")
            finally
                reporting = false
            end
        end
    end
end

# sliders send the design while they move: at most every 150 ms
@throttle design 150

# ---- client-side helpers (run in the browser) ---------------------------------------------

@methods """
part: function() { return this.nav === 'temp' ? this.tvar : this.nav },
spec: function(sect) { return (this.specs[this.nav] || {})[sect] || [] },
getv: function(sect, key) { const d = this.design[this.part()]; return d ? d[sect][key] : 0 },
// every edit replaces the whole design, so the change always reaches the server (an
// element changed inside one of its arrays did not)
edit: function(f) { const d = JSON.parse(JSON.stringify(this.design)); f(d[this.part()]); this.design = d },
setv: function(sect, key, v) {
    const x = Number(v);
    if (v === '' || v === null || isNaN(x)) return;
    this.edit(d => { d[sect][key] = x });
},
setPath: function(v) { this.edit(d => { d.target.path = v }) },
getm: function(key, m) { const d = this.design[this.part()]; return d ? d.month[key][m] : 0 },
setm: function(key, m, v) {
    const x = Number(v);
    if (v === '' || v === null || isNaN(x)) return;
    this.edit(d => { d.month[key][m] = x });
},
monthChanged: function(m) { return this.spec('month').some(c => this.getm(c.key, m) !== 0) },
// a month back to the values the design started with (0 for typical values, the file's departures otherwise)
startm: function(key, m) { const d = this.design0[this.part()]; return d ? d.month[key][m] : 0 },
monthEdited: function(m) { return this.spec('month').some(c => this.getm(c.key, m) !== this.startm(c.key, m)) },
resetMonth: function(m) { this.edit(d => this.spec('month').forEach(c => { d.month[c.key][m] = this.startm(c.key, m) })) },
resetMonths: function() { this.edit(d => Object.keys(d.month).forEach(k => { d.month[k] = d.month[k].map(() => 0) })) },
resetTargets: function() { this.edit(d => { Object.keys(d.target).forEach(k => { if (k !== 'path') d.target[k] = 0 }); Object.keys(d.mdelta).forEach(k => { d.mdelta[k] = d.mdelta[k].map(() => 0) }) }) },
targetsSet: function() { const t = this.design[this.part()].target; return Object.keys(t).some(k => k !== 'path' && t[k] !== 0) || this.mdeltaSet() },
// the monthly Δ: an extra change of each month by the end year, on top of the whole-year
// change, for the statistic picked above the table (md_key)
mdc: function() { const s = this.spec('mdelta'); return s.find(c => c.key === this.md_key) || s[0] || {} },
mdOpts: function() { return this.spec('mdelta').map(c => ({label: c.label + (this.mdeltaSet(c.key) ? ' •' : ''), value: c.key})) },
getd: function(m) { const d = this.design[this.part()], k = this.mdc().key; return d && d.mdelta && d.mdelta[k] ? d.mdelta[k][m] : 0 },
setd: function(m, v) {
    const x = Number(v), k = this.mdc().key;
    if (v === '' || v === null || isNaN(x)) return;
    this.edit(d => { d.mdelta[k][m] = x });
},
// any monthly Δ set (of statistic k, or of any)
mdeltaSet: function(k) {
    const d = this.design[this.part()];
    if (!d || !d.mdelta) return false;
    return (k ? [k] : Object.keys(d.mdelta)).some(j => (d.mdelta[j] || []).some(x => x !== 0));
},
resetMdelta: function() { this.edit(d => Object.keys(d.mdelta).forEach(k => { d.mdelta[k] = d.mdelta[k].map(() => 0) })) },
// the month's whole change by the end year: whole-year change plus its monthly Δ
totalDelta: function(m) {
    const c = this.mdc(), t = this.getv('target', c.tkey), x = this.getd(m);
    if (c.key === 'dr') { const p = ((1 + t / 100) * (1 + x / 100) - 1) * 100; return (p > 0 ? '+' : '') + Math.round(p) + '%' }
    const s = +(t + x).toFixed(3); return (s > 0 ? '+' : '') + s;
},
// ---- goals on the End year tab ----
// days above hot_thr: the goal typed, or the end year as designed
hotGoal: function() { return this.hot_goal >= 0 ? this.hot_goal : this.hot_days.end },
// the statistic whose trend chart is shown: the one picked, or the first in the table
trendStat: function() { return this.trend_charts[this.trend_stat] ? this.trend_stat : (this.trend_rows[0] || {}).stat },
setGoalThr: function(k, v) { const x = Number(v); if (v !== '' && v !== null && !isNaN(x)) this[k] = x },
setHotThr: function(v) { const x = Number(v); if (v !== '' && v !== null && !isNaN(x)) this.hot_thr = x },
setHotGoal: function(v) { const x = Number(v); if (v !== '' && v !== null && !isNaN(x) && x >= 0) this.hot_goal = x },
// the rain pattern of the whole year: start year, end year as designed, and the goal
rainNow: function(end) {
    const w = this.design.rain.whole, t = this.design.rain.target;
    const pww = w.pww + (end ? t.pww : 0), pwd = w.pwd + (end ? t.pwd : 0), total = w.total * (end ? 1 + t.total / 100 : 1);
    return this.rainPattern(pww, pwd, total);
},
rainPattern: function(pww, pwd, total) {
    // the chances kept within 0.01–0.99, as the designer keeps them (R_LIMITS), so the spell
    // lengths stay finite and positive whatever the change
    const lim = v => Math.min(0.99, Math.max(0.01, v));
    pww = lim(pww); pwd = lim(pwd);
    const wf = pwd / (1 - pww + pwd);
    return {pww: pww, pwd: pwd, total: total, dry: 1 / pwd, wet: 1 / (1 - pww), days: 365 * wf, mm: total / (365 * wf)};
},
rsDry: function() { return this.rs_dry > 0 ? this.rs_dry : +this.rainNow(true).dry.toFixed(1) },
rsMm: function() { return this.rs_mm > 0 ? this.rs_mm : +this.rainNow(true).mm.toFixed(1) },
setRs: function(k, v) { const x = Number(v); if (v !== '' && v !== null && !isNaN(x) && x > 0) this[k] = x },
clearRs: function() { this.rs_dry = -1; this.rs_mm = -1 },
rsSet: function() { return this.rs_dry > 0 || this.rs_mm > 0 },
// pwd from the dry spell (1 / pwd), then pww from the share of rainy days the rain per
// rainy day needs (share = pwd / (1 - pww + pwd)); null with the reason if it cannot be
rainGoal: function() {
    const e = this.rainNow(true), pwd = 1 / this.rsDry(), wf = e.total / (365 * this.rsMm());
    if (wf > 1) return {err: 'That rain per rainy day is too small for the annual rainfall: it would need rain every day.'};
    const pww = 1 + pwd - pwd / wf, w = this.design.rain.whole;
    if (pwd < 0.01 || pwd > 0.99) return {err: 'That dry spell is outside what the model allows (1 to 100 days).'};
    if (pww < 0.01 || pww > 0.99) return {err: 'These two goals cannot go together: ' + (pww < 0.01 ? 'rain would have to fall on too few days after dry spells this long; shorten the dry spell or lower the rain per rainy day.' : 'wet spells would never end; lengthen the dry spell or raise the rain per rainy day.')};
    const dww = +(pww - w.pww).toFixed(2), dwd = +(pwd - w.pwd).toFixed(2);
    if (Math.abs(dww) > 0.5 || Math.abs(dwd) > 0.5) return {err: 'That needs pww or pwd to change by more than 0.5, the most the End year allows.'};
    return Object.assign(this.rainPattern(w.pww + dww, w.pwd + dwd, e.total), {dww: dww, dwd: dwd});
},
applyRain: function() {
    const g = this.rainGoal();
    if (g.err) return;
    this.edit(d => { d.target.pww = g.dww; d.target.pwd = g.dwd });
    this.rs_dry = -1; this.rs_mm = -1;
},
// the end-year value a target change gives, as shown beside it
endValue: function(c) {
    const w = this.getv('whole', c.key), t = this.getv('target', c.key);
    if (c.key === 'peak' || c.key === 'dry') return 'day ' + (((Math.round(w + t) - 1) % 365 + 365) % 365 + 1);
    if (this.nav === 'rain' && c.key === 'total') return Math.round(w * (1 + t / 100)) + ' mm';
    return +(w + t).toFixed(3) + (c.unit && c.unit !== '% points' ? ' ' + c.unit : c.unit === '% points' ? ' %' : '');
},
dayLabel: function(d) {
    const dt = new Date(Date.UTC(2001, 0, Math.round(d)));
    return dt.toLocaleDateString('en-GB', {day: 'numeric', month: 'short', timeZone: 'UTC'});
},
years: function() { return Math.max(0, this.end_year - this.start_year + 1) },
siteSteps: function() { return [
    {nav: 'temp', no: 2, label: 'Temp', hint: 'Minimum and maximum air temperature'},
    {nav: 'wind', no: 3, label: 'Wind', hint: 'Wind speed'},
    {nav: 'rain', no: 4, label: 'Rain', hint: 'Rainfall and wet days'},
    {nav: 'stats', no: 5, label: 'Review', icon: 'fact_check', hint: 'Check the design, then save it as the stats file for MsiaGen'}] },
toolSteps: function() { return [
    {nav: 'run', icon: 'play_circle', label: 'Generate', hint: 'Run MsiaGen for the site, and check the fit'},
    {nav: 'simw', icon: 'search', label: 'Explore', hint: 'See what the generated weather looks like'}] },
showSiteSteps: function() { return this.nav !== 'home' && this.nav !== 'about' && !!this.site },
openSite: function(s) {
    if (this.design_dirty && !confirm('Open ' + s + '? The changes to the current design will be lost.')) return;
    this.opening = s;      // the card's button spins until the Site page opens
    this.open_site = s;
},
askDelete: function(s) { this.del_typed = ''; this.del_step = 1; this.del_target = s },
// step 1 to 2 only when the name is typed exactly (case sensitive)
nameTyped: function() { if (this.del_typed === this.del_target) this.del_step = 2 },
askRestart: function() { if (confirm('Start the design again from the chosen weather? Every change will be lost.')) this.restart_design = true },
askSave: function() {
    if (this.stats_exists && !confirm(this.statsName() + ' already exists in the site folder. Replace it?')) return;
    this.save_stats = true;
},
askRun: function() {
    const s = this.statsName();
    // asked only when there is a stats file to replace
    if (this.run_src === 'design' && this.run_has.stats && !confirm('Generating from the design replaces ' + s + ' with the design. Go on?')) return;
    if (this.run_src === 'obs' && this.run_has.stats && !confirm('Rebuilding from observed weather replaces ' + s + '. Go on?')) return;
    this.start_run = true;
},
statsName: function() { return (this.site || 'site').trim().replace(/\\s+/g, '-') + '-stats.csv' },
downloadStats: function() {
    const blob = new Blob([this.stats_text], {type: 'text/csv'});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = this.statsName(); a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
},
// the stats file as a table: latitude, header, rows, and the column groups (13 columns per
// statistic: the annual value, then the 12 months); parsed once per version of the file
statsTable: function() {
    if (this._st && this._st.text === this.stats_text) return this._st.t;
    const lines = (this.stats_text || '').split('\\n').filter(l => l.trim() !== '');
    const t = {lat: lines[0] || '', header: [], rows: [], groups: []};
    if (lines.length > 1) {
        t.header = lines[1].split(',');
        t.rows = lines.slice(2).map(l => l.split(','));
        const name = k => { const m = /^(mean|sd|rlag|skew)_(tmin|tmax|wind)\$/.exec(k);
            return m ? m[2] + ' · ' + m[1] : k === 'totrain' ? 'rain · total (mm)' : 'rain · ' + k };
        for (const h of t.header.slice(1)) {
            const k = h.replace(/\\d+\$/, '');
            const g = t.groups[t.groups.length - 1];
            if (g && g.key === k) g.n++; else t.groups.push({key: k, label: name(k), n: 1});
        }
    }
    this._st = {text: this.stats_text, t: t};
    return t;
},
baseUploaded: function(info) {
    const r = info.xhr.responseText;
    if (r.startsWith('ok')) this.upload_rev = this.upload_rev + 1;
    else this.toast = {msg: r, type: 'negative', t: Date.now()};
    this.\$refs.bup && this.\$refs.bup.reset();
},
figUrl: function(kind) { return '/figure?site=' + encodeURIComponent(this.fit_site) + '&kind=' + kind + '&r=' + this.figrev },
simUrl: function() { return '/download/sim?site=' + encodeURIComponent(this.statsName().replace(/-stats\\.csv\$/, '')) },
// the start- and end-year colours of the variable shown
yc: function() { return this.year_colors[this.part()] || ['#2f6aa3', '#d9661f'] },
gradeIcon: function(g) { return g === 'good' ? 'check_circle' : g === 'fair' ? 'check_circle_outline' : g === 'acceptable' ? 'error_outline' : 'cancel' },
plotConfig: function() { return {displaylogo: false, responsive: true, toImageButtonOptions: {format: 'png', scale: 2},
    showSendToCloud: false, showLink: false, showEditInChartStudio: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'sendDataToCloud', 'editInChartStudio']} }
"""

@watch """
nav: function() { this.opening = '' },
toast: function(v) {
    if (!v.msg) return;
    clearTimeout(this._toastTimer);
    this._toastTimer = setTimeout(() => { this.toast = {msg: '', type: v.type, t: 0} }, v.type === 'negative' ? 9000 : 4500);
}
"""

# Plotly throws when it resizes a chart that was removed in the meantime (charts are drawn
# afresh on every change); skip those, as Sawit does.
@mounted """
if (window.Plotly && !Plotly.Plots.resize.__msiagen) {
    const resize = Plotly.Plots.resize;
    Plotly.Plots.resize = function(gd) {
        const el = typeof gd === 'string' ? document.getElementById(gd) : gd;
        if (!el || !el.isConnected || el.offsetParent === null) return Promise.resolve();
        return resize.apply(this, arguments);
    };
    Plotly.Plots.resize.__msiagen = true;
}
"""

# ---- routes for files ----------------------------------------------------------------------

function routes()
    # Stipple's theme without its @import of a Google font (msiagen.css sets the fonts), so
    # the page works offline
    theme = Genie.Assets.asset_path(Stipple.assets_config, :css, path=Stipple.THEMES_FOLDER, file="theme-default-light")
    css = replace(read(joinpath(pkgdir(Stipple), "assets", "css", Stipple.THEMES_FOLDER, "theme-default-light.css"), String),
                  r"^@import url\([^)]*\);\s*"m => "")
    Genie.Router.route(theme) do
        Genie.Renderer.HTTP.Response(200, ["Content-Type" => "text/css"]; body=css)
    end
    Genie.Router.route("/download/sim") do
        s = Genie.Requests.getpayload(:site, "")
        (s in sites() && isfile(simfile(s))) || return Genie.Renderer.respond("No simulated weather.", 404)
        Genie.Renderer.HTTP.Response(200, ["Content-Type" => "text/csv",
            "Content-Disposition" => "attachment; filename=\"$(s)-sim.csv\""]; body=read(simfile(s), String))
    end
    Genie.Router.route("/figure") do
        s = Genie.Requests.getpayload(:site, "")
        k = Genie.Requests.getpayload(:kind, "")
        (s in sites() && k in ("sim", "fit") && isfile(figfile(s, k))) || return Genie.Renderer.respond("No figure.", 404)
        Genie.Renderer.HTTP.Response(200, ["Content-Type" => "image/png", "Cache-Control" => "no-store"]; body=read(figfile(s, k)))
    end
    Genie.Router.route("/upload/baseline", method=Genie.Router.POST) do
        files = Genie.Requests.filespayload()
        isempty(files) && return Genie.Renderer.respond("No file received.", 200)
        f = first(values(files))
        length(f.data) > MAX_UPLOAD && return Genie.Renderer.respond("$(f.name) is over $(MAX_UPLOAD ÷ 1_000_000) MB, too large for a daily weather file.", 200)
        text = String(copy(f.data))
        try
            average_year(parse_weather(text))
        catch e
            return Genie.Renderer.respond("Could not read $(f.name): $(sprint(showerror, e))", 200)
        end
        UPLOAD[] = (name=f.name, text=text)
        Genie.Renderer.respond("ok", 200)
    end
end

const LAYOUT = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <% Stipple.sesstoken() %>
    <title>MsiaGen</title>
    <link rel="icon" href="/favicon-32.png" sizes="32x32" type="image/png">
    <link rel="icon" href="/favicon.svg" type="image/svg+xml">
    <link rel="apple-touch-icon" href="/apple-touch-icon.png">
    <style>[v-cloak] { display: none; }</style>
    <% join(Stipple.Layout.stipplecore_css(), "\\n    ") %>
    <% join(Stipple.Layout.theme(; core_theme = false), "\\n    ") %>
    <link rel="stylesheet" href="/msiagen.css">
  </head>
  <body>
    <% Stipple.page(model, partial = true, v__cloak = true, [Stipple.Genie.Renderer.Html.@yield], Stipple.@if(:isready), core_theme = true, include_themes = false, class = "") %>
  </body>
</html>
"""

function port_in_use(port)
    try
        close(Sockets.listen(Sockets.localhost, port))
        false
    catch
        true
    end
end

"""
Run the page computations once on a site of the data folder before the page opens, so the
first site the user opens does not wait for Julia to compile them (about 5 s otherwise).
"""
function warmup()
    s = something(findfirst(x -> isfile(obsfile(x)) && isfile(statsfile(x)) && isfile(simfile(x)), sites()), 0)
    s == 0 && return
    site = sites()[s]
    b = average_year(parse_weather(read(obsfile(site), String)))
    d = new_design(b)
    d["tmin"]["target"]["M"] = 1.0
    st = all_stats(d, 2000:2002)
    for v ∈ TVARS
        x = st.t[v]
        overlay_chart([(x[1], t_sample_year(x[1], t_draws(0)), "#000000", "a"), (x[end], t_sample_year(x[end], t_draws(0)), "#000000", "b")];
                      unit="°C", ylabel="", other=[(x[1], "#000000", "c")])
    end
    overlay_chart([(st.w[1], w_sample_year(st.w[1], w_draws(0)), "#000000", "a")]; unit="m/s", ylabel="", floor_at=W_FLOOR)
    rain_overlay_chart([(st.r[1], "#000000", "a", 1), (st.r[end], "#000000", "b", 2)])
    rain_grid_chart(st.r[1], make_strip(st.r[1], r_draws(0)))
    value_checks(st); stats_csv(3.0, 2000:2002, st); trend_charts(st, 2000:2002)
    df = read_sim(site)
    sim_summary(df)
    for v ∈ ("tmin", "tmax", "wind", "rain")
        sim_page(site, v, 0, (hot=33.0, dry_day=0.5, dry_month=100.0))
    end
    fit_check(site)
    nothing
end

function start(; port=8100, open_browser=true)
    if port_in_use(port)
        println("\nMsiaGen is already running in another window (port $port).",
                "\nClose that window first, then start MsiaGen again.\n")
        exit(1)
    end
    println("  Preparing the computations...")
    try
        warmup()
    catch e
        @warn "warm-up skipped" exception = e
    end
    println("  Starting the web page...")
    routes()
    Genie.config.server_document_root = joinpath(dirname(@__DIR__), "ui", "public")
    template = joinpath(dirname(@__DIR__), "ui", "app.jl.html")
    # Genie renders only tags it knows; register every custom tag in the template up front
    for tag in unique(m.captures[1] for m in eachmatch(r"<([a-z]+(?:-[a-z]+)+)", read(template, String)))
        Genie.Renderer.Html.register_normal_element(replace(tag, "-" => "__"); context=Genie.Renderer.Html)
    end
    @page("/", template, layout=LAYOUT)
    up(port; async=true)
    url = "http://localhost:$port"
    println("  Preparing the screens (almost done)...")
    # (in a task of its own, so the server keeps answering while this waits)
    warm = @async try
        Genie.Renderer.HTTP.get(url; request_timeout=300, retry=false)
    catch
    end
    timedwait(() -> istaskdone(warm), 300.0; pollint=0.5)
    println("""

      ===========================================================
         MSIAGEN IS READY at $url
         Keep this window open while you use MsiaGen.
         To stop it, press Ctrl+C here or close this window.
      ===========================================================
    """)
    if open_browser
        Sys.isapple() ? run(`open $url`; wait=false) : Sys.iswindows() ? run(`cmd /c start $url`; wait=false) : nothing
    end
    url
end

end
