# Plotly charts (data and layout as Dicts, drawn by StipplePlotly in the browser).
# Colours: validated with the dataviz palette checks against the card surface #f9f8fb.

const INK, INK2, INK3 = "#1f1d26", "#4d4a57", "#625f6b"
const GRID, RULE = "#e7e5ec", "#d8d6de"
const PRIME, OCHRE = "#3a2d5c", "#8a5f17"
const C_TMIN, C_TMAX, C_WIND = "#2f6aa3", "#c0472f", "#1f8a70"
const T_COLOR = Dict("tmin" => C_TMIN, "tmax" => C_TMAX)
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

"""
Chart of one year's illustrative days `pts` around its monthly statistics `s` (μ, σ,
mdays, ndays): the ±1 sd band and mean of each month; optionally `ref`, another year's
monthly means (dashed), `floor_at`, a lower limit (dotted), and `other` = (s, label,
color), another variable's statistics drawn faded behind.
"""
function day_chart(s, pts; color, label, unit, ylabel, ref=nothing, floor_at=nothing,
                   floor_label="", other=nothing, height=360, yrange=nothing, ref_label="start-year monthly mean",
                   legend=true)
    starts = cumsum([1; s.mdays[1:end-1]])
    traces = Dict{String, Any}[]
    if !isnothing(other)
        o = other.s
        append!(traces, band_traces(starts, s.mdays, o.μ .- o.σ, o.μ .+ o.σ, other.color, "$(other.label) ±1 sd"; opacity=0.08))
        xo, yo = month_steps(starts, s.mdays, o.μ)
        push!(traces, Dict{String, Any}("x" => xo, "y" => yo, "mode" => "lines", "name" => "$(other.label) monthly mean",
            "line" => Dict("color" => rgba(other.color, 0.55), "width" => 2), "hoverinfo" => "skip"))
    end
    lo = isnothing(floor_at) ? s.μ .- s.σ : max.(floor_at, s.μ .- s.σ)
    append!(traces, band_traces(starts, s.mdays, lo, s.μ .+ s.σ, color, "$label ±1 sd"))
    xm, ym = month_steps(starts, s.mdays, s.μ)
    tips = repeat(["$(MONTH_NAMES[m]): mean $(r2(s.μ[m])) $unit, sd $(r2(s.σ[m])) $unit" for m ∈ 1:12]; inner=2)
    push!(traces, Dict{String, Any}("x" => xm, "y" => ym, "mode" => "lines", "name" => "$label monthly mean",
        "line" => Dict("color" => color, "width" => 2.5), "text" => tips, "hovertemplate" => "%{text}<extra></extra>"))
    if !isnothing(ref)
        xr, yr = month_steps(starts, s.mdays, ref)
        push!(traces, Dict{String, Any}("x" => xr, "y" => yr, "mode" => "lines", "name" => ref_label,
            "line" => Dict("color" => "#8d8a96", "width" => 1.6, "dash" => "dash"), "hoverinfo" => "skip"))
    end
    day_tips = ["$(MONTH_NAMES[m]) $k: $(r2(pts[starts[m] + k - 1])) $unit" for m ∈ 1:12 for k ∈ 1:s.mdays[m]]
    push!(traces, Dict{String, Any}("x" => collect(1:s.ndays), "y" => pts, "mode" => "lines+markers",
        "name" => "illustrative $label day", "text" => day_tips, "hovertemplate" => "%{text}<extra></extra>",
        "line" => Dict("color" => rgba(color, 0.35), "width" => 0.8),
        "marker" => Dict("color" => color, "size" => 5, "line" => Dict("color" => "#f9f8fb", "width" => 0.6))))
    shapes = month_lines(starts)
    isnothing(floor_at) || push!(shapes, Dict("type" => "line", "xref" => "paper", "x0" => 0, "x1" => 1,
        "y0" => floor_at, "y1" => floor_at, "line" => Dict("color" => INK3, "width" => 1.5, "dash" => "dot")))
    annots = isnothing(floor_at) ? [] : [Dict("xref" => "paper", "x" => 1, "y" => floor_at, "yanchor" => "bottom",
        "xanchor" => "right", "text" => floor_label, "showarrow" => false, "font" => Dict("size" => 12, "color" => INK3))]
    # a shared range, so the start- and end-year charts compare at a glance
    yr = !isnothing(yrange) ? Dict{String, Any}("range" => yrange) :
         isnothing(floor_at) ? Dict{String, Any}() : Dict{String, Any}("rangemode" => "tozero")
    Dict("data" => traces, "layout" => base_layout("height" => height, "shapes" => shapes, "annotations" => annots,
        "showlegend" => legend, "margin" => Dict("l" => 60, "r" => 16, "t" => 10, "b" => legend ? 48 : 30),
        "xaxis" => month_axis(starts, s.mdays, s.ndays),
        "yaxis" => axis("title" => Dict("text" => ylabel, "font" => Dict("size" => 13, "color" => INK2)), yr...)))
end

t_chart(s, pts, v; ref=nothing, other=nothing, kw...) =
    day_chart(s, pts; color=T_COLOR[v], label=v, unit="°C", ylabel="Air temperature (°C)", ref, kw...,
              other=isnothing(other) ? nothing : (s=other, label=other_tvar(v), color=T_COLOR[other_tvar(v)]))

w_chart(s, pts; ref=nothing, kw...) =
    day_chart(s, pts; color=C_WIND, label="wind", unit="m/s", ylabel="Wind speed (m/s)", ref, kw...,
              floor_at=W_FLOOR, floor_label="0.1 m/s floor")

"Illustrative daily rain (left, one row per month) and monthly rainfall (right) of year `s`."
function r_chart(s, strip; ref=nothing, top=nothing, ref_label="start-year total")
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
    scale = [[[(k - 1) / n, c], [k / n, c]] for (k, c) ∈ enumerate(RAIN_COLORS)] |> x -> reduce(vcat, x)
    traces = Dict{String, Any}[
        Dict("type" => "heatmap", "z" => z, "x" => collect(1:31), "y" => MONTH_NAMES, "text" => tips,
             "hovertemplate" => "%{text}<extra></extra>", "colorscale" => scale, "zmin" => -0.5, "zmax" => n - 0.5,
             "showscale" => false, "xgap" => 2, "ygap" => 2, "xaxis" => "x", "yaxis" => "y"),
        Dict("type" => "bar", "orientation" => "h", "x" => s.total, "y" => MONTH_NAMES, "xaxis" => "x2", "yaxis" => "y",
             "marker" => Dict("color" => RAIN_COLORS[5]), "name" => "monthly rainfall", "showlegend" => false,
             "hovertemplate" => "%{y}: %{x:.0f} mm<extra></extra>")]
    top = something(top, maximum(s.total))
    # each bar's label, placed clear of the bar and of the start-year tick
    edge = isnothing(ref) ? s.total : max.(s.total, ref)
    push!(traces, Dict("type" => "scatter", "mode" => "text", "x" => edge .+ 0.03 * max(1.0, top), "y" => MONTH_NAMES,
        "xaxis" => "x2", "yaxis" => "y", "text" => ["$(round(Int, s.total[m])) mm · $(s.raindays[m]) d" for m ∈ 1:12],
        "textposition" => "middle right", "textfont" => Dict("size" => 12, "color" => INK2), "cliponaxis" => false,
        "hoverinfo" => "skip", "showlegend" => false))
    if !isnothing(ref)
        top = max(top, maximum(ref))
        push!(traces, Dict("type" => "scatter", "mode" => "markers", "x" => ref, "y" => MONTH_NAMES, "xaxis" => "x2",
            "yaxis" => "y", "name" => ref_label, "marker" => Dict("symbol" => "line-ns-open", "size" => 20,
            "color" => "#8d8a96", "line" => Dict("width" => 2.5, "color" => "#8d8a96")), "hovertemplate" => "$ref_label: %{x:.0f} mm<extra></extra>"))
    end
    Dict("data" => traces, "layout" => base_layout("height" => 400, "bargap" => 0.3, "showlegend" => false,
        "margin" => Dict("l" => 44, "r" => 24, "t" => 28, "b" => 48),
        "xaxis" => axis("domain" => [0, 0.62], "tickvals" => [1; 5:5:30], "showgrid" => false, "ticks" => "",
                        "title" => Dict("text" => "Day of month", "font" => Dict("size" => 13, "color" => INK2))),
        "xaxis2" => axis("domain" => [0.68, 1], "range" => [0, 1.6 * max(1.0, top)], "anchor" => "y",
                         "title" => Dict("text" => "Rainfall (mm)", "font" => Dict("size" => 13, "color" => INK2))),
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

"Simulated daily weather of one year: temperatures, wind, rain (three stacked charts)."
function sim_year_charts(df, year)
    d = df[df.year .== year, :]
    x = [Date(year, m, dd) for (m, dd) ∈ zip(d.month, d.day)]
    xs = string.(x)
    lay(h, unit, kw::Pair...) = base_layout("height" => h, "margin" => Dict("l" => 56, "r" => 12, "t" => 8, "b" => 28),
        "xaxis" => axis("type" => "date", "tickformat" => "%b"), "yaxis" => axis("title" => Dict("text" => unit)),
        "hovermode" => "x unified", kw...)
    line(y, n, c) = Dict{String, Any}("x" => xs, "y" => y, "name" => n, "mode" => "lines",
                                     "line" => Dict("color" => c, "width" => 1.5), "hovertemplate" => "%{y:.1f}")
    charts = Any[]
    ("tmax" in names(d)) && push!(charts, Dict("title" => "Air temperature", "data" => [line(d.tmin, "tmin", C_TMIN), line(d.tmax, "tmax", C_TMAX)],
        "layout" => lay(260, "°C", "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.2))))
    ("wind" in names(d)) && push!(charts, Dict("title" => "Wind speed", "data" => [line(d.wind, "wind", C_WIND)],
        "layout" => lay(200, "m/s", "showlegend" => false)))
    ("rain" in names(d)) && push!(charts, Dict("title" => "Rainfall", "data" => [Dict{String, Any}("type" => "bar", "x" => xs, "y" => d.rain,
        "name" => "rain", "marker" => Dict("color" => RAIN_COLORS[5]), "hovertemplate" => "%{y:.1f} mm")],
        "layout" => lay(200, "mm", "showlegend" => false, "bargap" => 0.15)))
    charts
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

"The illustrative wet and dry days of one year: a grid, one row per month (no bars)."
function rain_grid_chart(s, strip; height=330)
    c = r_chart(s, strip)
    c["data"] = [first(c["data"])]
    c["layout"]["xaxis"]["domain"] = [0, 1]
    delete!(c["layout"], "xaxis2")
    c["layout"]["height"] = height
    c["layout"]["margin"] = Dict("l" => 44, "r" => 10, "t" => 6, "b" => 44)
    c
end
