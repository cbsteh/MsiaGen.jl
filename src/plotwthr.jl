# Charts of simulated daily weather for a site, all on one page (PNG), plus
# a per-year summary table (also printed). Built to stay readable for
# 1 to 30+ years: calendar-style charts give each year its own row of fixed
# height, so the page grows evenly with the number of years.
#
# Kept in its own submodule so plotting names (CairoMakie) cannot clash with
# names used by the weather generator (e.g. Distributions).

module Plotting

using CairoMakie, CSV, DataFrames, Dates, Printf, Statistics
import ..MONTHS

export plot_weather


# ---------------------------------------------------------------------------
# Colours (reference palette of the dataviz guidelines)
# ---------------------------------------------------------------------------

const SURFACE = colorant"#fcfcfb"
const INK = colorant"#0b0b0b"
const INK2 = colorant"#52514e"
const MUTED = colorant"#898781"
const GRID = colorant"#e1e0d9"
const NEUTRAL = colorant"#f0efec"
const BLUE = colorant"#2a78d6"          # categorical slot 1
const DARK_BLUE = colorant"#104281"     # sequential blue, step 650
const ORANGE = colorant"#eb6834"        # categorical slot 2

# rain classes: dry, then sequential blue (light = little, dark = a lot)
const RAIN_COLORS = [NEUTRAL, colorant"#cde2fb", colorant"#86b6ef",
                     colorant"#3987e5", colorant"#1c5cab", colorant"#0d366b"]

# wind classes: sequential orange (light = calm, dark = windy); orange
# because blue is already used for rain
const WIND_COLORS = [colorant"#fbe3d8", colorant"#f6bfa5", colorant"#f19670",
                     colorant"#eb6834", colorant"#b84c22", colorant"#7c3315"]
const WIND_BOUNDS = [1.0, 1.5, 2.0, 3.0, 4.0]      # m/s

# diverging, red end first: red = dry / hot, grey = at threshold, blue = wet / cool
const DIVERGING = [colorant"#b52f2f", colorant"#e34948", colorant"#f2a19c",
                   NEUTRAL,
                   colorant"#9ec5f4", colorant"#3987e5", colorant"#184f95"]

# Calendar charts place each day by its date in a leap year (`xday`), so
# every month starts at the same x in every year; 29 February stays empty
# in other years
const MONTH_START = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336]

const WIDTH = 1300      # width of a daily (calendar) plot area, px
const ROWH = 20         # height of one year in calendar charts, px
const CELLH = 24        # height of one year in the monthly grid, px
const MONTHLY_WIDTH = 12 * 56     # width of the monthly rain grid, px


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

function read_weather(fname)
    df = DataFrame(CSV.File(fname))
    df.date = Date.(df.year, df.month, df.day)
    df.x = xday.(df.date)
    sort!(df, :date)
end


# x position of a date in the calendar charts: its day of a leap year
xday(d::Date) = dayofyear(Date(2000, month(d), day(d)))


# All runs of consecutive dry days of at least `minlen` days, as (start, stop)
# dates. Missing days end a spell.
function find_dry_spells(df, threshold, minlen)
    spells = Tuple{Date,Date}[]
    start = nothing
    prev = nothing
    for (d, r) ∈ zip(df.date, df.rain)
        dry = !ismissing(r) && r < threshold
        gap = !isnothing(prev) && d - prev > Day(1)
        if !isnothing(start) && (!dry || gap)
            (prev - start).value + 1 >= minlen && push!(spells, (start, prev))
            start = nothing
        end
        dry && isnothing(start) && (start = d)
        prev = d
    end
    if !isnothing(start) && (prev - start).value + 1 >= minlen
        push!(spells, (start, prev))
    end
    spells
end


# Per-year indicators. A dry spell is counted in the year it starts.
function yearly_summary(df, spells, th)
    s = map(sort(unique(df.year))) do y
        d = df[df.year .== y, :]
        mrain = [sum(skipmissing(d.rain[d.month .== m]); init=0.0) for m ∈ 1:12]
        has_m = [any(d.month .== m) for m ∈ 1:12]
        lens = [(b - a).value + 1 for (a, b) ∈ spells if year(a) == y]
        (year=y,
         days=nrow(d),
         rain_mm=round(sum(skipmissing(d.rain)); digits=1),
         rain_days=count(r -> !ismissing(r) && r >= th.dry_day, d.rain),
         dry_months=count(mrain[has_m] .< th.dry_month),
         dry_spells=length(lens),
         longest_spell_days=isempty(lens) ? 0 : maximum(lens),
         hot_days=hasproperty(d, :tmax) ? count(t -> !ismissing(t) && t >= th.hot_day, d.tmax) : missing,
         mean_tmax=col_mean(d, :tmax, 1),
         mean_tmin=col_mean(d, :tmin, 1),
         mean_wind=col_mean(d, :wind, 2))
    end
    DataFrame(s)
end


# Mean of column `col` of `d`, rounded; missing if there is no such column
col_mean(d, col, digits) = hasproperty(d, col) ? round(mean(skipmissing(d[!, col])); digits=digits) : missing


# 0-based class of `x`, given ascending lower class bounds (first is -Inf)
classify(x, bounds) = ismissing(x) ? NaN : Float64(searchsortedlast(bounds, x) - 1)


fmt(x) = @sprintf("%g", x)


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

# A titled section of the page: title and subtitle on top; the chart goes
# in gl[3, 1] and its legend (if any) in gl[4, 1]
function section(pos, title, subtitle; valign=:center)
    gl = GridLayout(pos; valign=valign)
    Label(gl[1, 1], title; fontsize=16, font=:bold, color=INK, halign=:left,
          tellwidth=false)
    Label(gl[2, 1], subtitle; color=INK2, halign=:left, tellwidth=false)
    rowgap!(gl, 1, 2)
    gl
end


function style_axis!(ax)
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.leftspinecolor = GRID
    ax.bottomspinecolor = GRID
    ax.topspinevisible = false
    ax.rightspinevisible = false
    ax.xtickcolor = GRID
    ax.ytickcolor = GRID
    ax.xticklabelcolor = INK2
    ax.yticklabelcolor = INK2
    ax.xlabelcolor = INK2
    ax.ylabelcolor = INK2
    ax
end


function month_axis!(ax)
    ax.xticks = (MONTH_START .+ 14, MONTHS)
    xlims!(ax, 0.5, 366.5)
    vlines!(ax, MONTH_START[2:end] .- 0.5; color=GRID, linewidth=1)
end


# Axis with one row per year (earliest at the top)
function year_axis(pos, years, width, rowh)
    ny = length(years)
    ax = Axis(pos; width=width, height=rowh * ny, yreversed=true,
              yticks=(1:ny, string.(years)))
    ylims!(ax, ny + 0.5, 0.5)
    style_axis!(ax)
end


# One-row legend placed under a chart, so its size does not depend on how
# many years (rows) the chart has
function class_legend(pos, colors, labels, title)
    Legend(pos, [PolyElement(color=c, strokecolor=GRID, strokewidth=0.5) for c ∈ colors],
           labels, title; orientation=:horizontal, titleposition=:left,
           framevisible=false, labelcolor=INK2, titlecolor=INK2, halign=:left,
           tellwidth=false, colgap=14, padding=(0, 0, 0, 0))
end


# Daily values as a year × day grid, coloured by class
function calendar!(ax, df, years, col, bounds, colors)
    grid = fill(NaN, 366, length(years))
    yi = Dict(y => i for (i, y) ∈ enumerate(years))
    for (x, y, v) ∈ zip(df.x, df.year, df[!, col])
        grid[x, yi[y]] = classify(v, bounds)
    end
    n = length(colors)
    heatmap!(ax, 1:366, 1:length(years), grid;
             colormap=cgrad(colors, n; categorical=true),
             colorrange=(-0.5, n - 0.5), nan_color=:transparent)
    month_axis!(ax)
    length(years) > 1 && hlines!(ax, 1.5:1:length(years)-0.5; color=SURFACE, linewidth=2)
end


# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------

function rain_calendar!(pos, df, years, th)
    gl = section(pos, "Daily rain",
                 "Each row is a year, each cell a day. Grey: dry day " *
                 "(rain < $(fmt(th.dry_day)) mm).")
    ax = year_axis(gl[3, 1], years, WIDTH, ROWH)
    bounds = [th.dry_day, 5, 10, 20, 50]
    calendar!(ax, df, years, :rain, [-Inf; bounds], RAIN_COLORS)
    class_legend(gl[4, 1], RAIN_COLORS,
                 ["dry", "$(fmt(th.dry_day))–5", "5–10", "10–20", "20–50", "≥ 50"],
                 "mm/day")
end


function tmax_calendar!(pos, df, years, th)
    h = th.hot_day
    gl = section(pos, "Daily maximum temperature",
                 "Each row is a year, each cell a day. Red: hot day " *
                 "(tmax ≥ $(fmt(h)) °C); blue: cooler.")
    ax = year_axis(gl[3, 1], years, WIDTH, ROWH)
    bounds = [h - 3, h - 2, h - 1, h, h + 1, h + 2]
    calendar!(ax, df, years, :tmax, [-Inf; bounds], reverse(DIVERGING))
    class_legend(gl[4, 1], reverse(DIVERGING),
                 ["< $(fmt(h - 3))", "$(fmt(h - 3))–$(fmt(h - 2))",
                  "$(fmt(h - 2))–$(fmt(h - 1))", "$(fmt(h - 1))–$(fmt(h))",
                  "$(fmt(h))–$(fmt(h + 1))", "$(fmt(h + 1))–$(fmt(h + 2))",
                  "≥ $(fmt(h + 2))"],
                 "tmax (°C)")
end


function wind_calendar!(pos, df, years)
    gl = section(pos, "Daily wind speed",
                 "Each row is a year, each cell a day. Darker: windier.")
    ax = year_axis(gl[3, 1], years, WIDTH, ROWH)
    calendar!(ax, df, years, :wind, [-Inf; WIND_BOUNDS], WIND_COLORS)
    b = fmt.(WIND_BOUNDS)
    class_legend(gl[4, 1], WIND_COLORS,
                 ["< $(b[1])"; ["$(b[i])–$(b[i+1])" for i ∈ 1:length(b)-1]; "≥ $(b[end])"],
                 "m/s")
end


function dry_spell_chart!(pos, years, spells, th)
    gl = section(pos, "Dry spells of $(th.dry_spell) days or more",
                 "Each bar is one spell of consecutive dry days " *
                 "(rain < $(fmt(th.dry_day)) mm); label gives its length in days.")
    ax = year_axis(gl[3, 1], years, WIDTH, ROWH)
    month_axis!(ax)
    bars = Rect2f[]
    labels = []
    for (a, b) ∈ spells
        n = (b - a).value + 1
        pieces = []
        d = a
        while d <= b         # draw a spell crossing New Year on both rows
            e = min(b, Date(year(d), 12, 31))
            yi = findfirst(==(year(d)), years)
            x0, x1 = xday(d) - 0.5, xday(e) + 0.5
            push!(bars, Rect2f(x0, yi - 0.35, x1 - x0, 0.7))
            push!(pieces, (x0, x1, yi, n))
            d = e + Day(1)
        end
        push!(labels, argmax(p -> p[2] - p[1], pieces))   # label the longer part
    end
    isempty(bars) || poly!(ax, bars; color=ORANGE)
    # labels last, so no bar covers them: inside the bar if it is long
    # enough, otherwise beside it (on the left near the end of the year).
    # One text plot per alignment.
    placed = map(labels) do (x0, x1, yi, n)
        x, ha = x1 - x0 >= 10 ? ((x0 + x1) / 2, :center) :
                x1 > 350 ? (x0 - 2, :right) : (x1 + 2, :left)
        (Point2f(x, yi), string(n), ha)
    end
    for ha ∈ (:center, :right, :left)
        sel = filter(p -> p[3] == ha, placed)
        isempty(sel) || text!(ax, first.(sel); text=[p[2] for p ∈ sel], align=(ha, :center),
                              color=INK, fontsize=11)
    end
    isempty(spells) && text!(ax, 183, (length(years) + 1) / 2; text="no dry spells",
                             align=(:center, :center), color=MUTED)
end


function monthly_rain_chart!(pos, df, years, th)
    ny = length(years)
    gl = section(pos, "Monthly rain (mm)",
                 "Red: dry month (below $(fmt(th.dry_month)) mm); " *
                 "blue: above. Darker = further from the threshold."; valign=:top)
    mm = fill(NaN, 12, ny)
    yi = Dict(y => j for (j, y) ∈ enumerate(years))
    for g ∈ groupby(df, [:year, :month])
        mm[g.month[1], yi[g.year[1]]] = sum(skipmissing(g.rain))
    end
    ax = year_axis(gl[3, 1], years, MONTHLY_WIDTH, CELLH)
    ax.xticks = (1:12, MONTHS)
    xlims!(ax, 0.5, 12.5)
    lr = clamp.(log2.(max.(mm, 1.0) ./ th.dry_month), -2, 2)   # ×¼ … ×4 of threshold
    lr[isnan.(mm)] .= NaN
    heatmap!(ax, 1:12, 1:ny, lr; colormap=cgrad(DIVERGING), colorrange=(-2, 2),
             nan_color=:transparent)
    cells = [(m, j) for j ∈ 1:ny, m ∈ 1:12 if !isnan(mm[m, j])]
    text!(ax, [Point2f(m, j) for (m, j) ∈ cells];
          text=[string(round(Int, mm[m, j])) for (m, j) ∈ cells],
          color=[abs(lr[m, j]) > 1.2 ? colorant"#ffffff" : INK for (m, j) ∈ cells],
          align=(:center, :center), fontsize=11)
    ny > 1 && hlines!(ax, 1.5:1:ny-0.5; color=SURFACE, linewidth=2)
    vlines!(ax, 1.5:1:11.5; color=SURFACE, linewidth=2)
end


function cumulative_rain_chart!(pos, df, years)
    gl = section(pos, "Cumulative rain by year (mm)",
                 "Thin lines: each year; thick line: median across years."; valign=:top)
    height = max(CELLH * length(years), 300)     # at least as tall as the monthly grid
    ax = style_axis!(Axis(gl[3, 1]; width=WIDTH - MONTHLY_WIDTH - 90, height=height))
    ax.xticks = (MONTH_START[1:2:end] .+ 14, MONTHS[1:2:end])
    month_axis!(ax)
    ax.ygridvisible = true
    ax.ygridcolor = GRID
    # each full year's rain to date at every x (1-366); a year without
    # 29 February keeps its 28 February value there
    curves = Vector{Float64}[]
    for y ∈ years
        d = df[df.year .== y, :]
        c = cumsum(coalesce.(d.rain, 0.0))
        lines!(ax, d.x, c; color=(BLUE, length(years) > 10 ? 0.3 : 0.6),
               linewidth=1.2)
        if nrow(d) >= 365
            full = fill(NaN, 366)
            full[d.x] .= c
            isnan(full[60]) && (full[60] = full[59])
            push!(curves, full)
        end
    end
    if length(curves) >= 2
        med = [median(c[i] for c ∈ curves) for i ∈ 1:366]
        lines!(ax, 1:366, med; color=DARK_BLUE, linewidth=3)
        axislegend(ax, [LineElement(color=BLUE, linewidth=1.2),
                        LineElement(color=DARK_BLUE, linewidth=3)],
                   ["each year", "median"]; position=:lt, framevisible=false,
                   labelcolor=INK2)
    end
end


# Per-year summary table, as a grid of text labels
function summary_table!(pos, summ)
    gl = section(pos, "Summary by year", "")
    tbl = gl[3, 1] = GridLayout(; halign=:left, tellwidth=false)
    hdr = ["Year", "Rain (mm)", "Rain days", "Dry months", "Dry spells",
           "Longest spell (d)", "Hot days", "Mean tmax (°C)", "Mean tmin (°C)",
           "Mean wind (m/s)"]
    for (k, h) ∈ enumerate(hdr)
        Label(tbl[1, k], h; color=INK2, font=:bold, halign=:right)
    end
    for (i, row) ∈ enumerate(summary_rows(summ))
        for (k, v) ∈ enumerate(row)
            Label(tbl[i+1, k], v; color=INK, halign=:right)
        end
    end
    colgap!(tbl, 28)
    rowgap!(tbl, 3)
end


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

"""
    plot_weather(site; folder="data", dry_day=1.0, dry_spell=14,
                 dry_month=100.0, hot_day=33.0, show=true)

Chart the simulated daily weather in `folder/site/<site>-sim.csv` on one page,
saved as `<site>-sim.png` in the same folder:

- daily rain, one row per year
- dry spells of at least `dry_spell` days
- monthly rain against the dry-month threshold, and rain accumulated through
  each year
- daily maximum temperature, one row per year
- daily wind speed, one row per year
- a per-year summary table, also printed

The figure is displayed if `show=true`, and returned.

Thresholds:
- `dry_day`: a day with rain below this (mm) is a dry day
- `dry_spell`: dry spells of at least this many days are shown
- `dry_month`: a month with rain below this (mm) is a dry month
- `hot_day`: a day with tmax at or above this (°C) is a hot day
"""
function plot_weather(
    site::AbstractString;
    folder::AbstractString="data",
    dry_day::Real=1.0,
    dry_spell::Int=14,
    dry_month::Real=100.0,
    hot_day::Real=33.0,
    show::Bool=true
)
    path = joinpath(folder, site)
    df = read_weather(joinpath(path, "$(site)-sim.csv"))
    th = (; dry_day, dry_spell, dry_month, hot_day)
    years = sort(unique(df.year))
    spells = find_dry_spells(df, dry_day, dry_spell)
    summ = yearly_summary(df, spells, th)

    fig = Figure(backgroundcolor=SURFACE, fontsize=13,
                 figure_padding=(24, 48, 16, 24))     # left, right, bottom, top
    Label(fig[1, 1], "$(site): simulated daily weather, $(first(years))–$(last(years))";
          fontsize=22, font=:bold, color=INK, halign=:left, tellwidth=false)
    Label(fig[2, 1], "Dry day: rain < $(fmt(dry_day)) mm   ·   dry spell: ≥ $(dry_spell) " *
                     "dry days   ·   dry month: < $(fmt(dry_month)) mm   ·   " *
                     "hot day: tmax ≥ $(fmt(hot_day)) °C";
          color=INK2, halign=:left, tellwidth=false)

    rain_calendar!(fig[3, 1], df, years, th)
    dry_spell_chart!(fig[4, 1], years, spells, th)
    rain_row = GridLayout(fig[5, 1]; halign=:left)
    monthly_rain_chart!(rain_row[1, 1], df, years, th)
    cumulative_rain_chart!(rain_row[1, 2], df, years)
    colgap!(rain_row, 40)
    next = 6
    if hasproperty(df, :tmax)
        tmax_calendar!(fig[next, 1], df, years, th)
        next += 1
    end
    if hasproperty(df, :wind)
        wind_calendar!(fig[next, 1], df, years)
        next += 1
    end
    summary_table!(fig[next, 1], summ)
    rowgap!(fig.layout, 1, 4)
    rowgap!(fig.layout, 28)
    resize_to_layout!(fig)

    fname = joinpath(path, "$(site)-sim.png")
    save(fname, fig; px_per_unit=2)
    show && display(fig)

    println("Chart saved as $(fname)\n")
    show_summary(summ)
    fig
end


# `f(v)` as text, or a dash for a missing value
dash(f, v) = ismissing(v) ? "–" : f(v)


# One row of text per year for the summary table
summary_rows(summ) = [[string(s.year), @sprintf("%.0f", s.rain_mm),
                       string(s.rain_days), string(s.dry_months), string(s.dry_spells),
                       string(s.longest_spell_days), dash(string, s.hot_days),
                       dash(v -> @sprintf("%.1f", v), s.mean_tmax),
                       dash(v -> @sprintf("%.1f", v), s.mean_tmin),
                       dash(v -> @sprintf("%.2f", v), s.mean_wind)]
                      for s ∈ eachrow(summ)]


function show_summary(summ)
    hdr = ["Year", "Rain (mm)", "Rain days", "Dry months", "Dry spells",
           "Longest spell (d)", "Hot days", "Mean tmax", "Mean tmin", "Mean wind"]
    w = length.(hdr) .+ 2
    println(join(lpad.(hdr, w)))
    for row ∈ summary_rows(summ)
        println(join(lpad.(row, w)))
    end
end


# Specified against simulated statistics (from check_fit): one panel per
# statistic, one dot per month and year, with the 1:1 line
function plot_fit(site, pairs, fname; show::Bool=true)
    labels = unique(pairs.label)
    ncol = 4
    nrow_ = cld(length(labels), ncol)
    fig = Figure(backgroundcolor=SURFACE, fontsize=13,
                 figure_padding=(24, 32, 16, 24))
    Label(fig[0, 1:ncol], "$(site): simulated against specified statistics";
          fontsize=20, font=:bold, color=INK, halign=:left, tellwidth=false)
    Label(fig[1, 1:ncol], "One dot per month and year (annual total: per year). " *
                          "Dots on the dashed 1:1 line are a perfect fit. " *
                          "Rain pww and pwd come in steps, as a month's values are fractions of whole days.";
          color=INK2, halign=:left, tellwidth=false)
    for (i, lab) ∈ enumerate(labels)
        d = pairs[pairs.label .== lab, :]
        u = isempty(d.unit[1]) ? "" : " ($(d.unit[1]))"
        ax = style_axis!(Axis(fig[2 + (i - 1) ÷ ncol, 1 + (i - 1) % ncol];
                              width=230, height=200, title=lab, titlealign=:left,
                              titlesize=14, xlabel="specified" * u, ylabel="simulated" * u,
                              xgridvisible=true, ygridvisible=true,
                              xgridcolor=GRID, ygridcolor=GRID))
        # axis range: the data, but at least the statistic's minimum range
        lo, hi = extrema([d.spec; d.sim])
        mid, half = (lo + hi) / 2, max(hi - lo, maximum(d.min_span)) / 2
        lo, hi = mid - half, mid + half
        pad = max(1e-3, 0.05 * (hi - lo))
        lines!(ax, [lo - pad, hi + pad], [lo - pad, hi + pad];
               color=MUTED, linestyle=:dash, linewidth=1)
        scatter!(ax, d.spec, d.sim; color=(BLUE, 0.6), markersize=6)
        limits!(ax, lo - pad, hi + pad, lo - pad, hi + pad)
    end
    rowgap!(fig.layout, 0, 4)
    resize_to_layout!(fig)
    save(fname, fig; px_per_unit=2)
    show && display(fig)
    fig
end


end # module Plotting
