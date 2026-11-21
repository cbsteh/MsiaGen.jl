# Plotly charts (data and layout as Dicts, drawn by StipplePlotly in the browser).
# Colours: chosen to stay distinct with colour-blind vision, on the card surface #f9f8fb.

const INK, INK2, INK3 = "#1f1d26", "#4d4a57", "#625f6b"
const GRID, RULE = "#e7e5ec", "#d8d6de"
const PRIME, OCHRE = "#3a2d5c", "#8a5f17"
const C_TMIN, C_TMAX, C_WIND = "#2f6aa3", "#c0472f", "#1f8a70"
# daily rain classes: dry (neutral), then one blue hue from light to dark
const RAIN_COLORS = ["#e7e5ec", "#93b3db", "#6b94c9", "#4677b5", "#2b5a94", "#163b68"]
const RAIN_LABELS = ["dry", "< 5", "5–10", "10–20", "20–50", "≥ 50"]
const FONT = "Source Sans 3, -apple-system, Segoe UI, Roboto, Arial, sans-serif"

rgba(hex, a) = (h = lstrip(hex, '#'); "rgba($(parse(Int, h[1:2]; base=16)),$(parse(Int, h[3:4]; base=16)),$(parse(Int, h[5:6]; base=16)),$a)")
r2(x) = round(x; digits=2)

axis(kw::Pair...) = Dict{String, Any}("gridcolor" => GRID, "linecolor" => RULE, "zeroline" => false,
                                  "tickfont" => Dict("size" => 12.5, "color" => INK2),
                                  "title" => Dict("font" => Dict("size" => 13, "color" => INK2)), kw...)

function base_layout(kw::Pair...)
    Dict{String, Any}("height" => 420, "font" => Dict("family" => FONT, "size" => 13, "color" => INK),
        # the card's colour, not transparent: Plotly picks the toolbar's colours from it
        "paper_bgcolor" => "#f9f8fb", "plot_bgcolor" => "#f9f8fb",
        "modebar" => Dict("bgcolor" => "rgba(249,248,251,0.9)", "color" => INK3, "activecolor" => PRIME),
        "margin" => Dict("l" => 60, "r" => 16, "t" => 16, "b" => 48),
        "legend" => Dict("orientation" => "h", "x" => 0, "y" => -0.16, "font" => Dict("size" => 12.5, "color" => INK2)),
        "hoverlabel" => Dict("bgcolor" => "#ffffff", "bordercolor" => RULE, "font" => Dict("family" => FONT, "color" => INK)),
        "hovermode" => "closest", kw...)
end

# x and y of a step line over the months: flat across each month
function month_steps(starts, mdays, v)
    x, y = Float64[], Float64[]
    for m ∈ 1:12
        push!(x, starts[m] - 0.5, starts[m] + mdays[m] - 0.5)
        push!(y, v[m], v[m])
    end
    x, y
end

# a ±1 sd band (two traces: lower edge, then the upper edge filled down to it)
function band_traces(starts, mdays, lo, hi, color, name; opacity=0.14, legend=true)
    xl, yl = month_steps(starts, mdays, lo)
    xh, yh = month_steps(starts, mdays, hi)
    [Dict{String, Any}("x" => xl, "y" => yl, "mode" => "lines", "line" => Dict("width" => 0, "color" => color),
                       "hoverinfo" => "skip", "showlegend" => false),
     Dict{String, Any}("x" => xh, "y" => yh, "mode" => "lines", "line" => Dict("width" => 0, "color" => color),
                       "fill" => "tonexty", "fillcolor" => rgba(color, opacity), "hoverinfo" => "skip",
                       "name" => name, "showlegend" => legend)]
end

# month names under the middle of each month, thin lines between months
function month_axis(starts, mdays, ndays)
    axis("range" => [0.5, ndays + 0.5], "tickvals" => starts .+ (mdays .- 1) ./ 2,
         "ticktext" => MONTH_NAMES, "showgrid" => false, "ticks" => "")
end
month_lines(starts) = [Dict("type" => "line", "xref" => "x", "yref" => "paper", "x0" => s - 0.5, "x1" => s - 0.5,
                            "y0" => 0, "y1" => 1, "line" => Dict("color" => GRID, "width" => 1)) for s ∈ starts[2:end]]

"The illustrative wet and dry days of one year `s`: a grid, one row per month."
function rain_grid_chart(s, strip; height=330)
    z = Vector{Any}[]; tips = Vector{String}[]
    for m ∈ 1:12
        row, tip = Any[], String[]
        for d ∈ 1:31
            if d > s.mdays[m]
                push!(row, nothing); push!(tip, "")
            else
                r = strip[m][d]
                push!(row, r > 0 ? searchsortedlast(RAIN_BOUNDS, r) : 0)
                push!(tip, "$(MONTH_NAMES[m]) $d: " * (r > 0 ? "$(round(r; digits=1)) mm" : "dry"))
            end
        end
        push!(z, row); push!(tips, tip)
    end
    n = length(RAIN_COLORS)
    scale = reduce(vcat, [[[(k - 1) / n, c], [k / n, c]] for (k, c) ∈ enumerate(RAIN_COLORS)])
    trace = Dict("type" => "heatmap", "z" => z, "x" => collect(1:31), "y" => MONTH_NAMES, "text" => tips,
                 "hovertemplate" => "%{text}<extra></extra>", "colorscale" => scale, "zmin" => -0.5, "zmax" => n - 0.5,
                 "showscale" => false, "xgap" => 2, "ygap" => 2)
    Dict("data" => [trace], "layout" => base_layout("height" => height, "showlegend" => false,
        "margin" => Dict("l" => 44, "r" => 10, "t" => 6, "b" => 44),
        "xaxis" => axis("tickvals" => [1; 5:5:30], "showgrid" => false, "ticks" => "",
                        "title" => Dict("text" => "Day of month", "font" => Dict("size" => 13, "color" => INK2))),
        "yaxis" => axis("autorange" => "reversed", "showgrid" => false, "ticks" => "")))
end

"Small multiples of the annual values over the years: the path of the design."
function trend_charts(st, years)
    yrs = collect(years)
    mk(title, unit, series) = Dict("title" => title, "data" => [Dict{String, Any}("x" => yrs, "y" => y, "name" => n,
        "mode" => length(yrs) > 40 ? "lines" : "lines+markers", "line" => Dict("color" => c, "width" => 2),
        "marker" => Dict("size" => 6, "color" => c), "hovertemplate" => "%{x}: %{y:.2f} $unit<extra>$n</extra>")
        for (n, y, c) ∈ series],
        "layout" => base_layout("height" => 240, "showlegend" => length(series) > 1, "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.18),
            "margin" => Dict("l" => 52, "r" => 12, "t" => 10, "b" => 32),
            "xaxis" => axis("tickformat" => "d", "dtick" => max(1, ceil(Int, length(yrs) / 10))),
            "yaxis" => axis("title" => Dict("text" => unit), "range" => yrange(series))))
    annual(s) = sum(s.μ .* s.mdays) / s.ndays
    # a flat series gets a sensible axis, not one spread over its rounding error
    function yrange(series)
        v = reduce(vcat, [y for (_, y, _) ∈ series])
        lo, hi = extrema(v)
        pad = max(0.08 * (hi - lo), 0.02 * max(abs(hi), 1.0))
        [lo - pad, hi + pad]
    end
    [mk("Air temperature, annual mean", "°C", [("tmin", annual.(st.t["tmin"]), C_TMIN), ("tmax", annual.(st.t["tmax"]), C_TMAX)]),
     mk("Wind speed, annual mean", "m/s", [("wind", annual.(st.w), C_WIND)]),
     mk("Rainfall, annual total", "mm", [("rain", [sum(s.total) for s ∈ st.r], RAIN_COLORS[5])]),
     mk("Rain days per year", "days", [("rain days", [sum(s.raindays) for s ∈ st.r], RAIN_COLORS[4])])]
end


"""
Specified (stats file) against simulated values of one statistic: a dot per month and
year, the 1:1 line and the generator's fit tolerance around it. Dots outside it are ochre.
"""
function fit_chart(g, label, unit, tolf, min_span=0.0)
    spec, sim = Float64.(g.spec), Float64.(g.sim)
    tol = tolf.(spec)
    out = abs.(sim .- spec) .> tol
    # axis range: the data, but at least `min_span` (10 times the fit tolerance, as in
    # MsiaGen's own fit chart), so a small error does not fill the panel
    lo, hi = extrema([spec; sim])
    mid, half = (lo + hi) / 2, max(hi - lo, min_span, 1e-6) / 2
    lo, hi = mid - 1.06half, mid + 1.06half
    tips = [(g.month[i] == 0 ? "" : MONTH_NAMES[g.month[i]] * " ") * string(g.year[i]) *
            ": specified $(round(spec[i]; sigdigits=4)), simulated $(round(sim[i]; sigdigits=4))" for i ∈ eachindex(spec)]
    xs = range(lo, hi; length=60)
    band = [Dict{String, Any}("x" => collect(xs), "y" => collect(xs) .- tolf.(xs), "mode" => "lines", "hoverinfo" => "skip",
                              "line" => Dict("width" => 0), "showlegend" => false),
            Dict{String, Any}("x" => collect(xs), "y" => collect(xs) .+ tolf.(xs), "mode" => "lines", "hoverinfo" => "skip",
                              "line" => Dict("width" => 0), "fill" => "tonexty", "fillcolor" => rgba(PRIME, 0.08), "name" => "fit tolerance")]
    pts(sel, c, n) = Dict{String, Any}("x" => spec[sel], "y" => sim[sel], "mode" => "markers", "name" => n, "text" => tips[sel],
        "hovertemplate" => "%{text}<extra></extra>", "marker" => Dict("color" => c, "size" => 7, "opacity" => 0.85,
        "line" => Dict("color" => "#f9f8fb", "width" => 1)))
    u = isempty(unit) ? "" : " ($unit)"
    Dict("title" => label, "data" => [band; Dict{String, Any}("x" => [lo, hi], "y" => [lo, hi], "mode" => "lines", "name" => "1:1",
            "hoverinfo" => "skip", "line" => Dict("color" => INK3, "width" => 1.2, "dash" => "dash"));
            pts(.!out, PRIME, "within tolerance"); pts(out, OCHRE, "outside tolerance")],
         "layout" => base_layout("height" => 300, "showlegend" => false, "margin" => Dict("l" => 56, "r" => 12, "t" => 8, "b" => 44),
            "xaxis" => axis("range" => [lo, hi], "title" => Dict("text" => "specified" * u)),
            "yaxis" => axis("range" => [lo, hi], "title" => Dict("text" => "simulated" * u), "scaleanchor" => "x")))
end

# ---- start and end year on one chart -----------------------------------------------------

# the two years: blue and orange, a pair that stays distinct with colour-blind vision
const C_START, C_END = "#2f6aa3", "#d9661f"
# the start- and end-year colours of each variable: tmax has its own pair (purple, gold), so
# its lines never look like tmin's
const YEAR_COLORS = Dict("tmin" => (C_START, C_END), "tmax" => ("#6b4ea8", "#c99a1a"),
                         "wind" => (C_START, C_END), "rain" => (C_START, C_END))

"""
One chart of the start and end year together. `years`: (s, pts, color, label) of each year
shown (s: monthly statistics, pts: illustrative days); `other`: (s, color, label) of the
other temperature variable's years, drawn as dotted means behind; `days`: draw the
illustrative days.
"""
function overlay_chart(years; unit, ylabel, floor_at=nothing, floor_label="", other=[], days=true, height=440)
    traces = Dict{String, Any}[]
    vals = Any[]
    for (s, col, lab) ∈ other
        starts = cumsum([1; s.mdays[1:end-1]])
        # its ±1 sd band (lighter than the main variable's), then its mean, dotted
        for t ∈ band_traces(starts, s.mdays, s.μ .- s.σ, s.μ .+ s.σ, col, "$lab ±1 sd"; opacity=0.09)
            t["legendgroup"] = lab
            push!(traces, t)
        end
        x, y = month_steps(starts, s.mdays, s.μ)
        mt = repeat(["$lab, $(MONTH_NAMES[m]): mean $(r2(s.μ[m])) $unit, sd $(r2(s.σ[m])) $unit" for m ∈ 1:12]; inner=2)
        push!(traces, Dict{String, Any}("x" => x, "y" => y, "mode" => "lines", "name" => "$lab mean", "legendgroup" => lab,
            "line" => Dict("color" => col, "width" => 1.8, "dash" => "dot"), "opacity" => 0.8,
            "text" => mt, "hovertemplate" => "%{text}<extra></extra>"))
        append!(vals, [s.μ .- s.σ, s.μ .+ s.σ])
    end
    for (s, pts, col, lab) ∈ years
        starts = cumsum([1; s.mdays[1:end-1]])
        lo = isnothing(floor_at) ? s.μ .- s.σ : max.(floor_at, s.μ .- s.σ)
        for t ∈ band_traces(starts, s.mdays, lo, s.μ .+ s.σ, col, "$lab: ±1 sd"; opacity=0.13)
            t["legendgroup"] = lab
            push!(traces, t)
        end
        if days
            tips = ["$lab, $(MONTH_NAMES[m]) $k: $(r2(pts[starts[m] + k - 1])) $unit" for m ∈ 1:12 for k ∈ 1:s.mdays[m]]
            # joined day to day: how smooth or jumpy the line is shows rlag
            push!(traces, Dict{String, Any}("x" => collect(1:s.ndays), "y" => pts, "mode" => "lines+markers", "name" => "$lab: illustrative days",
                "legendgroup" => lab, "text" => tips, "hovertemplate" => "%{text}<extra></extra>",
                "line" => Dict("color" => rgba(col, 0.45), "width" => 1),
                "marker" => Dict("color" => rgba(col, 0.7), "size" => 4.5)))
            push!(vals, pts)
        end
        x, y = month_steps(starts, s.mdays, s.μ)
        mt = repeat(["$lab, $(MONTH_NAMES[m]): mean $(r2(s.μ[m])) $unit, sd $(r2(s.σ[m])) $unit" for m ∈ 1:12]; inner=2)
        push!(traces, Dict{String, Any}("x" => x, "y" => y, "mode" => "lines", "name" => "$lab: monthly mean", "legendgroup" => lab,
            "line" => Dict("color" => col, "width" => 3), "text" => mt, "hovertemplate" => "%{text}<extra></extra>"))
        append!(vals, [lo, s.μ .+ s.σ])
    end
    s0 = isempty(years) ? first(other)[1] : first(first(years))
    starts = cumsum([1; s0.mdays[1:end-1]])
    shapes = month_lines(starts)
    isnothing(floor_at) || push!(shapes, Dict("type" => "line", "xref" => "paper", "x0" => 0, "x1" => 1,
        "y0" => floor_at, "y1" => floor_at, "line" => Dict("color" => INK3, "width" => 1.5, "dash" => "dot")))
    annots = isnothing(floor_at) ? [] : [Dict("xref" => "paper", "x" => 1, "y" => floor_at, "yanchor" => "bottom",
        "xanchor" => "right", "text" => floor_label, "showarrow" => false, "font" => Dict("size" => 12, "color" => INK3))]
    yr = isnothing(floor_at) ? Dict{String, Any}() : Dict{String, Any}("rangemode" => "tozero")
    Dict("data" => traces, "layout" => base_layout("height" => height, "shapes" => shapes, "annotations" => annots,
        "legend" => Dict("orientation" => "h", "x" => 0, "y" => -0.13, "font" => Dict("size" => 12.5, "color" => INK2), "traceorder" => "grouped"),
        "margin" => Dict("l" => 60, "r" => 16, "t" => 10, "b" => 60),
        "xaxis" => month_axis(starts, s0.mdays, s0.ndays),
        "yaxis" => axis("title" => Dict("text" => ylabel, "font" => Dict("size" => 13, "color" => INK2)), yr...)))
end

# rain days are drawn as dots on stems, in their own colours, so they read apart from the
# rainfall bars
const RAINDAY_COLORS = ("#13867a", "#b8264f")    # start year (teal), end year (crimson)

"""
Monthly rainfall (grouped bars) and rain days (dots on stems) of the start and end year.
`years`: (s, color, label, which) with which = 1 for the start year, 2 for the end year.
"""
function rain_overlay_chart(years; height=440)
    traces = Dict{String, Any}[]
    n = length(years)
    for (s, col, lab) ∈ years
        push!(traces, Dict{String, Any}("type" => "bar", "orientation" => "h", "x" => s.total, "y" => MONTH_NAMES,
            "xaxis" => "x", "yaxis" => "y", "name" => "$lab: rainfall", "marker" => Dict("color" => col),
            "hovertemplate" => "$lab, %{y}: %{x:.0f} mm<extra></extra>"))
    end
    for (k, (s, _, lab, which)) ∈ enumerate(years)
        col = RAINDAY_COLORS[which]
        y = (1:12) .+ (n == 1 ? 0.0 : k == 1 ? -0.18 : 0.18)
        sx, sy = Any[], Any[]
        for m ∈ 1:12
            push!(sx, 0, s.raindays[m], nothing); push!(sy, y[m], y[m], nothing)
        end
        push!(traces, Dict{String, Any}("x" => sx, "y" => sy, "mode" => "lines", "xaxis" => "x2", "yaxis" => "y2",
            "line" => Dict("color" => col, "width" => 2), "hoverinfo" => "skip", "showlegend" => false))
        push!(traces, Dict{String, Any}("x" => s.raindays, "y" => collect(y), "mode" => "markers", "xaxis" => "x2", "yaxis" => "y2",
            "name" => "$lab: rain days", "marker" => Dict("color" => col, "size" => 9, "line" => Dict("color" => "#f9f8fb", "width" => 1)),
            "text" => MONTH_NAMES, "hovertemplate" => "$lab, %{text}: %{x} rain days<extra></extra>"))
    end
    title(t) = Dict("text" => t, "font" => Dict("size" => 13, "color" => INK2))
    Dict("data" => traces, "layout" => base_layout("height" => height, "barmode" => "group", "bargap" => 0.25, "bargroupgap" => 0.05,
        "legend" => Dict("orientation" => "h", "x" => 0, "y" => -0.13, "font" => Dict("size" => 12.5, "color" => INK2)),
        "margin" => Dict("l" => 44, "r" => 16, "t" => 10, "b" => 60),
        "xaxis" => axis("domain" => [0, 0.56], "rangemode" => "tozero", "title" => title("Monthly rainfall (mm)")),
        "xaxis2" => axis("domain" => [0.66, 1], "anchor" => "y2", "rangemode" => "tozero", "title" => title("Rain days")),
        "yaxis" => axis("autorange" => "reversed", "showgrid" => false, "ticks" => ""),
        # the month names again beside the rain days
        "yaxis2" => axis("anchor" => "x2", "range" => [12.6, 0.4], "tickvals" => collect(1:12), "ticktext" => MONTH_NAMES,
                         "showgrid" => false, "ticks" => "")))
end

# ---- the trend and goal checks ------------------------------------------------------------

"""
The trend check of one statistic: for each month, its change over the years in the stats
file (designed) and in the generated weather, with the band of chance (95%)
around the generated change; months off are drawn in ochre.
"""
function trend_chart(per, label, unit, y0, y1)
    u = isempty(unit) ? "" : " ($unit)"
    off = [p.off for p ∈ per]
    gen(sel, c, n) = Dict{String, Any}("x" => MONTH_NAMES[sel], "y" => [p.generated for p ∈ per[sel]], "mode" => "markers",
        "name" => n, "marker" => Dict("color" => c, "size" => 9, "symbol" => "circle"),
        "error_y" => Dict("type" => "data", "array" => [p.band for p ∈ per[sel]], "color" => rgba(c, 0.5), "thickness" => 1.5, "width" => 4),
        "text" => ["generated $(signs(p.generated)), trend error $(signs(p.err)) (chance ±$(fmt3(p.band)))" for p ∈ per[sel]],
        "hovertemplate" => "%{x}: %{text}<extra></extra>")
    Dict("title" => label,
         "data" => [Dict{String, Any}("x" => MONTH_NAMES, "y" => [p.designed for p ∈ per], "mode" => "markers", "name" => "designed",
                        "marker" => Dict("color" => INK3, "size" => 11, "symbol" => "line-ew-open", "line" => Dict("width" => 2.5, "color" => INK3)),
                        "hovertemplate" => "%{x}: designed %{y:.3g}<extra></extra>"),
                    gen(findall(.!off), PRIME, "generated"), gen(findall(off), OCHRE, "generated, off")],
         "layout" => base_layout("height" => 320, "margin" => Dict("l" => 60, "r" => 12, "t" => 8, "b" => 70),
            "xaxis" => axis("title" => Dict("text" => "")),
            "yaxis" => axis("title" => Dict("text" => "change $(y0)–$(y1)" * u), "zeroline" => true, "zerolinecolor" => RULE)))
end

"The goal check of one measure: expected (line) and generated (dots) in each year."
# `tol`: how far the averages may differ and still count as close. The y axis spans at least
# that much either side of the expected average, so small year-to-year differences do not
# fill the chart and look like a poor fit.
function goal_chart(years, e, g, label, unit, tol)
    v = filter(isfinite, [e; g])
    ef = filter(isfinite, e)
    c = isempty(ef) ? (isempty(v) ? 0.0 : sum(v) / length(v)) : sum(ef) / length(ef)
    lo, hi = min(minimum(v; init=c), c - tol), max(maximum(v; init=c), c + tol)
    pad = 0.05 * (hi - lo)
    # counts and lengths are never below 0; the axis still starts a pad below the lowest
    # value (below 0 too), so a dot at 0 is drawn whole, not cut in half by the axis
    lo = max(lo, 0.0)
    Dict("title" => label,
         "data" => [Dict{String, Any}("x" => years, "y" => e, "mode" => "lines", "name" => "expected",
                        "line" => Dict("color" => INK3, "width" => 2), "hovertemplate" => "%{x}: expected %{y:.3g}<extra></extra>"),
                    Dict{String, Any}("x" => years, "y" => g, "mode" => "markers", "name" => "generated",
                        "marker" => Dict("color" => PRIME, "size" => 7), "hovertemplate" => "%{x}: generated %{y:.3g}<extra></extra>")],
         "layout" => base_layout("height" => 330, "margin" => Dict("l" => 56, "r" => 12, "t" => 8, "b" => 70),
            "xaxis" => axis("title" => Dict("text" => ""), "tickformat" => "d"),
            "yaxis" => axis("title" => Dict("text" => unit), "range" => [lo - pad, hi + pad])))
end
