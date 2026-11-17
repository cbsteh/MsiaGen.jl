module WeatherApp

using Stipple, Stipple.ReactiveTools, StippleUI, StipplePlotly
using CSV, DataFrames, Dates, JSON, Printf, Random, Sockets, Statistics
using Distributions, SpecialFunctions
import MsiaGen
const Genie = Stipple.Genie
# Vue's <transition> (the step bar's sliding site group) is not an HTML element Genie knows;
# without this it warns UndefVarError(:transition) on the first page load
Genie.Renderer.Html.register_normal_element("transition"; context=Genie.Renderer.Html)

include("designer.jl")
include("charts.jl")
include("sites.jl")

# ---- the controls of each page (label, unit, range, step), drawn by the page template -------

ctl(key, label, unit, min, max, step; hint="", more=false) =
    Dict{String, Any}("key" => key, "label" => label, "unit" => unit, "min" => min, "max" => max,
                      "step" => step, "hint" => hint, "more" => more)

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
        "rcols" => ["rainfall (mm)", "pww", "pwd", "rain days", "wet spell (d)", "dry spell (d)"]))

const PATH_OPTS = [Dict("label" => l, "value" => v) for (v, l) ∈ PATHS]

# ---- the design and its baseline ------------------------------------------------------------

"The design of all four parts, from a file's average year `b` (or typical values: nothing)."
function new_design(b)
    get_b(k) = b isa NamedTuple ? getproperty(b, k) : nothing
    Dict{String, Any}("tmin" => temp_design("tmin", get_b(:tmin)), "tmax" => temp_design("tmax", get_b(:tmax)),
                      "wind" => wind_design(get_b(:wind)), "rain" => rain_design(get_b(:rain)))
end

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
    f = baseline_text(m.base_src[], m.site[])
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
    d = new_design(b)
    DESIGN0[] = d
    m.design[] = d
end
const DESIGN0 = Ref(new_design(nothing))   # the design as last started, to tell if it was changed

# ---- what each page shows -------------------------------------------------------------------

design_years(m) = m.start_year[]:max(m.start_year[], m.end_year[])
signed(v; digits=3) = (v = round(v; digits); v == 0 ? "" : (v > 0 ? "+" : "−") * string(abs(v)))
cell(v, d=0.0; digits=3) = Dict("v" => string(v), "d" => signed(d; digits))
capped_note(ys, part) = isempty(ys) ? String[] :
    ["Some $(part) values went past their limits in $(length(ys)) year(s), $(first(ys))–$(last(ys)). They were capped at the limit."]

function var_page!(m)
    nav = m.nav[]
    nav ∈ ("temp", "wind", "rain") || return
    st = all_stats(m.design[], design_years(m))
    pick(list) = m.view_end[] ? last(list) : first(list)
    if nav == "temp"
        v = m.tvar[]
        d = m.design[][v]["month"]
        s = pick(st.t[v])
        m.vchart[] = t_chart(s, t_sample_year(s, t_draws(m.t_new[])), v;
                             ref=m.view_end[] ? first(st.t[v]).μ : nothing, other=pick(st.t[other_tvar(v)]))
        m.vrows[] = [Dict("month" => MONTH_NAMES[k], "cells" => [cell(r2(s.μ[k]), d["dm"][k]), cell(s.σ[k], d["dsd"][k]),
                          cell(s.ρ[k], d["drl"][k]), cell(s.γ[k], d["dsk"][k])]) for k ∈ 1:12]
        m.vsummary[] = "Whole year $(s.year): annual mean $(r2(s.M)) °C, amplitude $(r2(s.A)) °C, warmest day $(show_day(s.peak))."
        n = t_cross_note(st.t)
        m.vnotes[] = [capped_note(capped_years(st).temperature, "temperature"); isnothing(n) ? String[] : [n]]
    elseif nav == "wind"
        d = m.design[]["wind"]["month"]
        s = pick(st.w)
        m.vchart[] = w_chart(s, w_sample_year(s, w_draws(m.w_new[])); ref=m.view_end[] ? first(st.w).μ : nothing)
        m.vrows[] = [Dict("month" => MONTH_NAMES[k], "cells" => [cell(r2(s.μ[k]), d["dm"][k]), cell(s.σ[k], d["dsd"][k]),
                          cell(s.ρ[k], d["drl"][k])]) for k ∈ 1:12]
        m.vsummary[] = "Whole year $(s.year): annual mean $(r2(s.M)) m/s, amplitude $(r2(s.A)) m/s, windiest day $(show_day(s.peak))."
        n = w_floor_note(st.w)
        m.vnotes[] = [capped_note(capped_years(st).wind, "wind"); isnothing(n) ? String[] : [n]]
    else
        d = m.design[]["rain"]["month"]
        s = pick(st.r)
        m.vchart[] = r_chart(s, make_strip(s, r_draws(m.r_new[])); ref=m.view_end[] ? first(st.r).total : nothing)
        m.vrows[] = [Dict("month" => MONTH_NAMES[k], "cells" => [cell(round(Int, s.total[k]), d["dr"][k]; digits=0),
                          cell(s.pww[k], d["dww"][k]), cell(s.pwd[k], d["dwd"][k]), cell(s.raindays[k]),
                          cell(round(s.wetspell[k]; digits=1)), cell(round(s.dryspell[k]; digits=1))]) for k ∈ 1:12]
        m.vsummary[] = "Whole year $(s.year): annual rainfall $(round(Int, sum(s.total))) mm, " *
                       "$(sum(s.raindays)) rain days, seasonality $(round(Int, s.seas))%, driest day $(show_day(s.dry))."
        m.vnotes[] = capped_note(capped_years(st).rain, "rainfall")
    end
    m.vrev[] += 1
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
    m.stats_exists[] = !isempty(strip(m.site[])) && isfile(statsfile(site_folder(m.site[])))
    m.srev[] += 1
end

"Problems that stop the stats file being written (empty if none)."
function design_problems(m)
    p = String[]
    isempty(strip(m.site[])) && push!(p, "Give the site a name (Site page).")
    !isempty(strip(m.site[])) && !valid_site(m.site[]) && push!(p, "The site name may hold letters, digits, spaces, '-', '_' and '.' only.")
    m.start_year[] > m.end_year[] && push!(p, "The start year is after the end year (Site page).")
    p
end

function save_stats!(m)
    s = site_folder(m.site[])
    mkpath(sitedir(s))
    years = design_years(m)
    write(statsfile(s), stats_csv(m.lat[], years, all_stats(m.design[], years)))
    m.site_cards[] = site_card.(sites())
    m.stats_exists[] = true
    s
end

function load_sim!(m, s)
    df = read_sim(s)
    if isnothing(df)
        m.sim_rows[] = Dict{String, Any}[]; m.sim_years[] = Int[]; m.sim_charts[] = Any[]
        return
    end
    m.sim_rows[] = sim_summary(df)
    ys = sort(unique(df.year))
    m.sim_years[] = ys
    m.sim_year[] in ys || (m.sim_year[] = first(ys))
    m.sim_charts[] = sim_year_charts(df, m.sim_year[])
    m.simrev[] += 1
end

function load_fit!(m, s)
    m.fit_site[] = s
    if !(isfile(statsfile(s)) && isfile(simfile(s)))
        m.fit_rows[] = Dict{String, Any}[]; m.fit_charts[] = Any[]
        m.fit_notes[] = ["$(s) needs both $(s)-stats.csv and $(s)-sim.csv. Generate its weather first."]
        return
    end
    try
        f = fit_check(s)
        m.fit_rows[] = f.rows; m.fit_charts[] = f.charts; m.fit_notes[] = f.notes
    catch e
        m.fit_rows[] = Dict{String, Any}[]; m.fit_charts[] = Any[]
        m.fit_notes[] = ["Could not check the fit: $(sprint(showerror, e))"]
    end
    m.figs[] = Dict("sim" => isfile(figfile(s, "sim")), "fit" => isfile(figfile(s, "fit")))
    m.fitrev[] += 1
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
    @in new_dialog = false
    @in new_name = ""
    @in do_new = false

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
    @in tvar = "tmin"
    @in ctab = "whole"
    @in view_end = false
    @in t_new = 0
    @in w_new = 0
    @in r_new = 0
    @out vchart = Dict("data" => [], "layout" => Dict())
    @out vrows = Dict{String, Any}[]
    @out vsummary = ""
    @out vnotes = String[]
    @out vrev = 0

    # statistics
    @out checks = String[]
    @out annual_rows = Dict{String, Any}[]
    @out stats_text = ""
    @out trends = Any[]
    @out stats_exists = false
    @out srev = 0
    @in save_stats = false

    # generate
    @in run_src = "design"
    @in seed = -1
    @in wetdry = MsiaGen.WETDRY_TMAX
    @in start_run = false
    @out running = false
    @out run_log = ""
    @out sim_rows = Dict{String, Any}[]
    @out sim_years = Int[]
    @in sim_year = 0
    @out sim_charts = Any[]
    @out simrev = 0

    # fit check
    @out fit_site = ""
    @out fit_rows = Dict{String, Any}[]
    @out fit_charts = Any[]
    @out fit_notes = String[]
    @out fitrev = 0
    @in make_reports = false
    @out reporting = false
    @out figs = Dict("sim" => false, "fit" => false)
    @out figrev = 0

    @onchange design begin
        design_dirty = design != DESIGN0[]
        var_page!(__model__)
        stats_page!(__model__)
    end
    @onchange nav, tvar, view_end, t_new, w_new, r_new, start_year, end_year begin
        var_page!(__model__)
        stats_page!(__model__)
        if nav == "fit" && !isempty(strip(site))
            s = site_folder(site)
            fit_site == s && !isempty(fit_rows) || load_fit!(__model__, s)
        end
    end
    @onchange lat, site begin
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
            base_opts = base_options(site)
            base_src == "upload" ? apply_baseline!(__model__) : (base_src = "upload")
        end
    end
    @onchange open_site begin
        if !isempty(open_site)
            open_site!(__model__, open_site)
            open_site = ""
            nav = "site"
        end
    end
    @onbutton do_new begin
        name = site_folder(new_name)
        if !valid_site(name)
            toast = toastmsg("A site name may hold letters, digits, spaces, '-', '_' and '.' only.", "negative")
        elseif name in sites()
            toast = toastmsg("$name already exists: open it from its card.", "warning")
        else
            site = name
            is_new = true
            base_opts = base_options(name)
            base_src == "typical" ? apply_baseline!(__model__) : (base_src = "typical")
            run_src = "design"
            load_sim!(__model__, name)
            new_dialog = false
            new_name = ""
            nav = "site"
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
        s = site_folder(site)
        p = run_src == "design" ? design_problems(__model__) : (isempty(strip(site)) ? ["Open or create a site first."] : String[])
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
                toast = toastmsg("Generated $(length(sim_years)) year(s) of daily weather for $s. Check the fit next.")
            catch e
                @error "generate" exception = (e, catch_backtrace())
                run_log = sprint(showerror, e)
                toast = toastmsg("MsiaGen stopped: see the messages below.", "negative")
            finally
                running = false
            end
        end
    end
    @onchange sim_year begin
        s = site_folder(site)
        if sim_year > 0 && isfile(simfile(s))
            sim_charts = sim_year_charts(read_sim(s), sim_year)
            simrev = simrev + 1
        end
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
resetMonth: function(m) { this.edit(d => this.spec('month').forEach(c => { d.month[c.key][m] = 0 })) },
resetMonths: function() { this.edit(d => Object.keys(d.month).forEach(k => { d.month[k] = d.month[k].map(() => 0) })) },
resetTargets: function() { this.edit(d => Object.keys(d.target).forEach(k => { if (k !== 'path') d.target[k] = 0 })) },
targetsSet: function() { const t = this.design[this.part()].target; return Object.keys(t).some(k => k !== 'path' && t[k] !== 0) },
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
fmtv: function(c, v) { return (c.key === 'peak' || c.key === 'dry') && c.unit === 'day of year' ? v + ' (' + this.dayLabel(v) + ')' : v },
signed: function(v) { return v > 0 ? '+' + v : String(v).replace('-', '−') },
years: function() { return Math.max(0, this.end_year - this.start_year + 1) },
siteSteps: function() { return [
    {nav: 'temp', no: 2, label: 'Temperature', hint: 'Minimum and maximum air temperature'},
    {nav: 'wind', no: 3, label: 'Wind', hint: 'Wind speed'},
    {nav: 'rain', no: 4, label: 'Rainfall', hint: 'Rainfall and wet days'},
    {nav: 'stats', no: 5, label: 'Statistics', hint: 'Checks, and the stats file for MsiaGen'}] },
toolSteps: function() { return [
    {nav: 'run', icon: 'play_circle', label: 'Generate', hint: 'Run MsiaGen for the site'},
    {nav: 'fit', icon: 'query_stats', label: 'Fit check', hint: 'Simulated against specified statistics'}] },
showSiteSteps: function() { return this.nav !== 'home' && this.nav !== 'about' && !!this.site },
openSite: function(s) {
    if (this.design_dirty && !confirm('Open ' + s + '? The changes to the current design will be lost.')) return;
    this.open_site = s;
},
askRestart: function() { if (confirm('Start the design again from the chosen weather? Every change will be lost.')) this.restart_design = true },
askSave: function() {
    if (this.stats_exists && !confirm(this.statsName() + ' already exists in the site folder. Replace it?')) return;
    this.save_stats = true;
},
askRun: function() {
    const s = this.statsName();
    if (this.run_src === 'design' && this.stats_exists && !confirm('Generating from the design replaces ' + s + ' with the design. Go on?')) return;
    if (this.run_src === 'obs' && !confirm('Rebuilding from observed weather replaces ' + s + '. Go on?')) return;
    this.start_run = true;
},
statsName: function() { return (this.site || 'site').trim().replace(/\\s+/g, '-') + '-stats.csv' },
downloadStats: function() {
    const blob = new Blob([this.stats_text], {type: 'text/csv'});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = this.statsName(); a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
},
csvPreview: function() {
    const lines = this.stats_text.split('\\n');
    return lines.slice(0, 14).join('\\n') + (lines.length > 15 ? '\\n…  (' + (lines.length - 3) + ' years in all)' : '');
},
baseUploaded: function(info) {
    const r = info.xhr.responseText;
    if (r.startsWith('ok')) this.upload_rev = this.upload_rev + 1;
    else this.toast = {msg: r, type: 'negative', t: Date.now()};
    this.\$refs.bup && this.\$refs.bup.reset();
},
figUrl: function(kind) { return '/figure?site=' + encodeURIComponent(this.fit_site) + '&kind=' + kind + '&r=' + this.figrev },
simUrl: function() { return '/download/sim?site=' + encodeURIComponent(this.statsName().replace(/-stats\\.csv\$/, '')) },
rainLegend: function() { return ['#e7e5ec', '#93b3db', '#6b94c9', '#4677b5', '#2b5a94', '#163b68'].map((c, i) =>
    ({color: c, label: ['dry', '< 5', '5–10', '10–20', '20–50', '≥ 50'][i]})) },
gradeIcon: function(g) { return g === 'good' ? 'check_circle' : g === 'fair' ? 'check_circle_outline' : g === 'acceptable' ? 'error_outline' : 'cancel' },
plotConfig: function() { return {displaylogo: false, responsive: true, toImageButtonOptions: {format: 'png', scale: 2},
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d']} }
"""

@watch """
// the chart shows the year the controls on the left are setting: the End year tab the end
// year, the Start year and Each month tabs the start year
ctab: function(v) { this.view_end = v === 'target' },
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
    <link rel="icon" href="/favicon.svg" type="image/svg+xml">
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

function start(; port=8100, open_browser=true)
    if port_in_use(port)
        println("\nMsiaGen is already running in another window (port $port).",
                "\nClose that window first, then start MsiaGen again.\n")
        exit(1)
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
