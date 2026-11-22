# The Explore tab of the Generate page: what the generated daily weather (<site>-sim.csv) looks like.
# Distributions, the monthly spread, day-to-day persistence, wet and dry spells, tmin
# against tmax, rain through the year, and the extremes of each year.

const SIM_CACHE = Ref((file="", mtime=0.0, df=DataFrame()))

"The simulated weather of site `s`, read again only when the file changes (nothing if none)."
function sim_weather(s)
    f = simfile(s)
    isfile(f) || return nothing
    t = mtime(f)
    c = SIM_CACHE[]
    (c.file == f && c.mtime == t) && return c.df
    df = sort!(DataFrame(CSV.File(f)), [:year, :month, :day])
    SIM_CACHE[] = (file=f, mtime=t, df=df)
    df
end

const SW_VARS = Dict("tmin" => (label="Min. air temp.", unit="°C", color=C_TMIN),
                     "tmax" => (label="Max. air temp.", unit="°C", color=C_TMAX),
                     "wind" => (label="Wind speed", unit="m/s", color=C_WIND),
                     "rain" => (label="Rainfall", unit="mm", color=RAIN_COLORS[5]))
const C_WET, C_DRY = "#2b5a94", "#b07d2b"

# ---- numbers ---------------------------------------------------------------------------------

"Gaussian kernel density of `x` at `n` points from its lowest to its highest value (Silverman's bandwidth)."
function kde(x; n=200)
    m = length(x)
    s = m > 1 ? std(x) : 0.0
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    h = 0.9 * max(min(s, iqr / 1.34), 1e-6) * m^(-0.2)
    # evaluated only over the data's own range: no curve beyond the lowest or highest value
    lo, hi = extrema(x)
    g = range(lo, hi; length=n)
    k = 1 / (m * h * sqrt(2π))
    collect(g), [k * sum(exp(-0.5 * ((gi - xi) / h)^2) for xi ∈ x) for gi ∈ g]
end

"Histogram of `x` as density: bin centres, densities and the bin width (Freedman–Diaconis bins)."
function density_hist(x)
    lo, hi = extrema(x)
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    w = 2 * max(iqr, 1e-9) * length(x)^(-1 / 3)
    nb = clamp(ceil(Int, (hi - lo) / max(w, 1e-9)), 10, 60)
    w = max(hi - lo, 1e-9) / nb
    counts = zeros(nb)
    for v ∈ x
        counts[clamp(floor(Int, (v - lo) / w) + 1, 1, nb)] += 1
    end
    lo .+ w .* ((1:nb) .- 0.5), counts ./ (length(x) * w), w
end

"Lengths of the runs of `true` in `flags`."
function runs(flags)
    out = Int[]
    n = 0
    for f ∈ flags
        if f
            n += 1
        elseif n > 0
            push!(out, n); n = 0
        end
    end
    n > 0 && push!(out, n)
    out
end

fmtn(x; d=2) = string(round(x; digits=d))

# ---- charts -----------------------------------------------------------------------------------

titled(t) = Dict("text" => t, "font" => Dict("size" => 13, "color" => INK2))

"Histogram and smoothed curve of `x`; `logx`: drawn on log10(x) with tick labels in x."
function dist_chart(x, color, xlabel; logx=false)
    v = logx ? log10.(x) : x
    c, d, w = density_hist(v)
    gx, gy = kde(v)
    lab(t) = logx ? "$(round(10^t; sigdigits=3))" : "$(round(t; sigdigits=4))"
    tips = ["$(lab(ci - w / 2))–$(lab(ci + w / 2)): $(round(100 * di * w; digits=1))% of days" for (ci, di) ∈ zip(c, d)]
    xa = logx ? axis("title" => titled(xlabel), "tickvals" => collect(-1:3), "ticktext" => ["0.1", "1", "10", "100", "1000"],
                     "range" => [minimum(gx), maximum(gx)]) : axis("title" => titled(xlabel))
    Dict("data" => [
            Dict("type" => "bar", "x" => c, "y" => d, "width" => w * 0.96, "marker" => Dict("color" => rgba(color, 0.35)),
                 "name" => "histogram", "text" => tips, "textposition" => "none", "hovertemplate" => "%{text}<extra></extra>"),
            Dict("x" => gx, "y" => gy, "mode" => "lines", "line" => Dict("color" => color, "width" => 2.5),
                 "name" => "smoothed (kernel density)", "hoverinfo" => "skip")],
         "layout" => base_layout("height" => 320, "bargap" => 0, "showlegend" => true,
            "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.12),
            "margin" => Dict("l" => 56, "r" => 12, "t" => 26, "b" => 48),
            "xaxis" => xa, "yaxis" => axis("title" => titled("density"), "rangemode" => "tozero")))
end

"""
One violin per month (the spread of its values), with a small box: middle half and median.
`logy`: drawn on log10 of the values (the violins' shapes then follow the log scale), with
tick labels in the values themselves.
"""
function month_violin_chart(groups, color, ylabel; logy=false)
    x = reduce(vcat, [fill(MONTH_NAMES[m], length(g)) for (m, g) ∈ enumerate(groups)])
    y = reduce(vcat, groups)
    logy && (y = log10.(y))
    # with few values per month (e.g. monthly totals of a few years) each value is also drawn
    few = maximum(length.(groups)) < 15
    tr = Dict{String, Any}("type" => "violin", "x" => x, "y" => y, "name" => ylabel,
        "points" => few ? "all" : false, "pointpos" => 0, "jitter" => 0.3,
        "marker" => Dict("color" => color, "size" => 5),
        "spanmode" => "hard", "scalemode" => "width", "width" => 0.85,
        "box" => Dict("visible" => true, "width" => 0.18, "fillcolor" => "#f9f8fb", "line" => Dict("color" => color, "width" => 1.2)),
        "meanline" => Dict("visible" => false),
        "line" => Dict("color" => color, "width" => 1.2), "fillcolor" => rgba(color, 0.25),
        "hoverinfo" => logy ? "skip" : "y", "hoveron" => "violins")
    # on the log scale the violin's own hover would give log values: a hidden point at each
    # month's median gives its numbers in the values themselves
    hov = Dict{String, Any}[]
    if logy
        ms = [m for m ∈ 1:12 if !isempty(groups[m])]
        q(m, p) = round(quantile(groups[m], p); sigdigits=3)
        push!(hov, Dict{String, Any}("type" => "scatter", "x" => MONTH_NAMES[ms], "y" => [log10(q(m, 0.5)) for m ∈ ms],
            "mode" => "markers", "marker" => Dict("size" => 14, "color" => "rgba(0,0,0,0)"), "showlegend" => false,
            "text" => ["$(MONTH_NAMES[m]): median $(q(m, 0.5)), middle half $(q(m, 0.25))–$(q(m, 0.75)), " *
                       "largest $(round(maximum(groups[m]); sigdigits=3)) ($(length(groups[m])) days)" for m ∈ ms],
            "hovertemplate" => "%{text}<extra></extra>"))
    end
    ticks = [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]
    ya = logy ? axis("title" => titled(ylabel), "tickvals" => log10.(ticks), "ticktext" => string.(ticks) .|> t -> replace(t, r"\.0$" => "")) :
                axis("title" => titled(ylabel))
    Dict("data" => [tr; hov], "layout" => base_layout("height" => 320, "showlegend" => false, "violinmode" => "overlay",
        "margin" => Dict("l" => 56, "r" => 12, "t" => 12, "b" => 36),
        "xaxis" => axis("showgrid" => false, "categoryorder" => "array", "categoryarray" => MONTH_NAMES),
        "yaxis" => ya))
end

"Each day against the day before, with the 1:1 line and the lag-1 correlation r."
function lag_chart(x, color, unit)
    a, b = x[1:end-1], x[2:end]
    # many days make a solid block: fainter dots the more there are
    op = length(a) > 5000 ? 0.12 : length(a) > 2000 ? 0.2 : 0.35
    lo, hi = extrema(x)
    pad = 0.04 * max(hi - lo, 1e-6)
    r = [lo - pad, hi + pad]
    Dict("data" => [
            Dict("type" => "scattergl", "x" => a, "y" => b, "mode" => "markers", "name" => "days",
                 "marker" => Dict("color" => rgba(color, op), "size" => 4), "hovertemplate" => "yesterday %{x:.2f}, today %{y:.2f} $unit<extra></extra>"),
            Dict("x" => [lo, hi], "y" => [lo, hi], "mode" => "lines", "name" => "1:1", "hoverinfo" => "skip",
                 "line" => Dict("color" => INK3, "width" => 1.2, "dash" => "dash"))],
         "layout" => base_layout("height" => 320, "showlegend" => false, "margin" => Dict("l" => 56, "r" => 12, "t" => 12, "b" => 48),
            "xaxis" => axis("title" => titled("yesterday ($unit)"), "range" => r),
            "yaxis" => axis("title" => titled("today ($unit)"), "range" => r),
            "annotations" => [Dict("x" => hi, "y" => hi, "text" => "1:1", "showarrow" => false, "xanchor" => "right", "yanchor" => "bottom",
                                   "xshift" => -4, "font" => Dict("size" => 12, "color" => INK3)),
                             Dict("xref" => "paper", "yref" => "paper", "x" => 0.02, "y" => 0.98, "xanchor" => "left", "yanchor" => "top",
                                  "text" => "r = " * replace(string(round(lag1(x); digits=2)), "-" => "−"), "showarrow" => false,
                                  "font" => Dict("size" => 13, "color" => INK))]))
end

"How often wet and dry spells of each length occur (spells longer than 30 days in one bar)."
function spell_chart(wet, dry)
    cnt(r) = [count(==(k), r) for k ∈ 1:30] |> v -> [v; count(>(30), r)]
    xs = [string.(1:30); ">30"]
    Dict("data" => [
            Dict("type" => "bar", "x" => xs, "y" => cnt(wet), "name" => "wet spells", "marker" => Dict("color" => C_WET),
                 "hovertemplate" => "%{y} wet spells of %{x} days<extra></extra>"),
            Dict("type" => "bar", "x" => xs, "y" => cnt(dry), "name" => "dry spells", "marker" => Dict("color" => C_DRY),
                 "hovertemplate" => "%{y} dry spells of %{x} days<extra></extra>")],
         "layout" => base_layout("height" => 320, "barmode" => "group", "bargap" => 0.2,
            "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.12), "margin" => Dict("l" => 56, "r" => 12, "t" => 26, "b" => 48),
            "xaxis" => axis("title" => titled("spell length (days)"), "type" => "category", "showgrid" => false),
            "yaxis" => axis("title" => titled("number of spells"))))
end

"""
pww and pwd of each month of each year: the share of wet days followed by a wet day, and of
dry days followed by a wet day, over the days whose next day is in the same month (as
MsiaGen counts them). Returns two 12 × years matrices (NaN where a month has no wet or no
dry day).
"""
function wet_chances(d, wet)
    ys = sort(unique(d.year))
    pww, pwd = fill(NaN, 12, length(ys)), fill(NaN, 12, length(ys))
    for (j, y) ∈ enumerate(ys), m ∈ 1:12
        w = wet[(d.year .== y) .& (d.month .== m)]
        length(w) < 2 && continue
        prev, next = w[1:end-1], w[2:end]
        nw, nd = count(prev), count(.!prev)
        nw > 0 && (pww[m, j] = count(prev .& next) / nw)
        nd > 0 && (pwd[m, j] = count(.!prev .& next) / nd)
    end
    pww, pwd
end

nanmean_row(A) = [(v = filter(!isnan, A[m, :]); isempty(v) ? NaN : mean(v)) for m ∈ 1:size(A, 1)]
nanext_row(A, f) = [(v = filter(!isnan, A[m, :]); isempty(v) ? NaN : f(v)) for m ∈ 1:size(A, 1)]

"pww and pwd through the year: each month's mean over the years, with a band for their range."
function chances_chart(pww, pwd)
    traces = Dict{String, Any}[]
    x = collect(1:12)
    for (A, col, name) ∈ ((pww, C_WET, "pww: wet day after a wet day"), (pwd, C_DRY, "pwd: wet day after a dry day"))
        mu = nanmean_row(A)
        if size(A, 2) > 1
            push!(traces, Dict("x" => x, "y" => nanext_row(A, minimum), "mode" => "lines", "line" => Dict("width" => 0, "color" => col),
                               "hoverinfo" => "skip", "showlegend" => false))
            push!(traces, Dict("x" => x, "y" => nanext_row(A, maximum), "mode" => "lines", "line" => Dict("width" => 0, "color" => col),
                               "fill" => "tonexty", "fillcolor" => rgba(col, 0.15), "hoverinfo" => "skip", "showlegend" => false))
        end
        push!(traces, Dict("x" => x, "y" => mu, "mode" => "lines+markers", "name" => name, "text" => MONTH_NAMES,
            "line" => Dict("color" => col, "width" => 2.5), "marker" => Dict("color" => col, "size" => 7),
            "hovertemplate" => "%{text}: %{y:.2f}<extra>$(split(name, ':')[1])</extra>"))
    end
    Dict("data" => traces, "layout" => base_layout("height" => 320,
        "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.12), "margin" => Dict("l" => 56, "r" => 12, "t" => 26, "b" => 36),
        "xaxis" => axis("tickvals" => x, "ticktext" => MONTH_NAMES, "showgrid" => false, "range" => [0.5, 12.5]),
        "yaxis" => axis("title" => titled("chance"), "range" => [0, 1])))
end

"Monthly mean tmax against monthly rainfall: one dot per month of each year."
function tmax_rain_chart(d)
    g = combine(groupby(d, [:year, :month]), :tmax => mean => :tmax, :rain => sum => :rain)
    tips = ["$(MONTH_NAMES[m]) $y: $(round(r; digits=0)) mm, tmax $(round(t; digits=2)) °C" for (y, m, r, t) ∈ zip(g.year, g.month, g.rain, g.tmax)]
    dots = Dict("type" => "scatter", "x" => g.rain, "y" => g.tmax, "mode" => "markers", "name" => "months",
            "text" => tips, "hovertemplate" => "%{text}<extra></extra>",
            "marker" => Dict("color" => rgba(C_TMAX, 0.6), "size" => 7, "line" => Dict("color" => "#f9f8fb", "width" => 1)))
    # least-squares line and its correlation r (none if rain or tmax does not vary)
    fit = nothing
    if nrow(g) >= 3 && std(g.rain) > 0 && std(g.tmax) > 0
        r = cor(g.rain, g.tmax)
        slope = r * std(g.tmax) / std(g.rain)
        x = [extrema(g.rain)...]
        y = mean(g.tmax) .+ slope .* (x .- mean(g.rain))
        fit = (line=Dict("x" => x, "y" => y, "mode" => "lines", "name" => "regression", "hoverinfo" => "skip",
                         "line" => Dict("color" => INK3, "width" => 1.5)),
               label=Dict("xref" => "paper", "yref" => "paper", "x" => 0.98, "y" => 0.96, "xanchor" => "right", "yanchor" => "top",
                          "text" => "r = " * replace(string(round(r; digits=2)), "-" => "−"), "showarrow" => false, "font" => Dict("size" => 13, "color" => INK)))
    end
    Dict("data" => isnothing(fit) ? [dots] : [dots, fit.line],
         "layout" => base_layout("height" => 320, "showlegend" => false, "margin" => Dict("l" => 56, "r" => 12, "t" => 12, "b" => 48),
            "xaxis" => axis("title" => titled("monthly rainfall (mm)"), "rangemode" => "tozero"),
            "yaxis" => axis("title" => titled("monthly mean tmax (°C)")),
            "annotations" => isnothing(fit) ? [] : [fit.label]))
end

"""
Rain accumulated through each year (day of year on x). With 5 years or more: each year a
faint line, with the median and the 10th to 90th percentile band of the years; the year
picked (`sel`, 0 for none) in bold. With fewer years: one line per year.
"""
function cumrain_chart(df, years, sel=0)
    traces = Dict{String, Any}[]
    n = length(years)
    starts = cumsum([1; [daysinmonth(2001, m) for m ∈ 1:11]])
    cum = Dict(y => cumsum(df.rain[df.year .== y]) for y ∈ years)
    if n < 5
        for (i, y) ∈ enumerate(years)
            col = n == 1 ? RAIN_COLORS[5] : RAIN_COLORS[2 + round(Int, 3 * (i - 1) / max(n - 1, 1))]
            push!(traces, Dict("x" => collect(1:length(cum[y])), "y" => cum[y], "mode" => "lines", "name" => string(y),
                "line" => Dict("color" => col, "width" => y == sel ? 3 : 1.6), "hovertemplate" => "$y, day %{x}: %{y:.0f} mm<extra></extra>"))
        end
    else
        for y ∈ years
            push!(traces, Dict("x" => collect(1:length(cum[y])), "y" => cum[y], "mode" => "lines", "name" => string(y), "showlegend" => false,
                "line" => Dict("color" => rgba(RAIN_COLORS[3], 0.3), "width" => 0.8), "hovertemplate" => "$y, day %{x}: %{y:.0f} mm<extra></extra>"))
        end
        # percentiles over the years, day by day (to day 365: a leap year's last day left out)
        M = reduce(hcat, [cum[y][1:365] for y ∈ years if length(cum[y]) >= 365])
        pct(p) = [quantile(view(M, k, :), p) for k ∈ 1:365]
        x = collect(1:365)
        push!(traces, Dict("x" => x, "y" => pct(0.1), "mode" => "lines", "line" => Dict("width" => 0), "hoverinfo" => "skip", "showlegend" => false))
        push!(traces, Dict("x" => x, "y" => pct(0.9), "mode" => "lines", "line" => Dict("width" => 0), "fill" => "tonexty",
            "fillcolor" => rgba(RAIN_COLORS[5], 0.18), "name" => "10th–90th percentile", "hoverinfo" => "skip"))
        push!(traces, Dict("x" => x, "y" => pct(0.5), "mode" => "lines", "name" => "median",
            "line" => Dict("color" => RAIN_COLORS[6], "width" => 2.5), "hovertemplate" => "median, day %{x}: %{y:.0f} mm<extra></extra>"))
        if sel ∈ years
            push!(traces, Dict("x" => collect(1:length(cum[sel])), "y" => cum[sel], "mode" => "lines", "name" => string(sel),
                "line" => Dict("color" => C_END, "width" => 2.5), "hovertemplate" => "$sel, day %{x}: %{y:.0f} mm<extra></extra>"))
        end
    end
    Dict("data" => traces, "layout" => base_layout("height" => 320, "showlegend" => n < 5 ? n <= 12 : true,
        "legend" => n < 5 ? Dict("orientation" => "v", "x" => 1.02, "y" => 1) : Dict("orientation" => "h", "x" => 0, "y" => 1.12),
        "margin" => Dict("l" => 60, "r" => 12, "t" => n < 5 ? 12 : 26, "b" => 36),
        "xaxis" => axis("tickvals" => starts .+ 14, "ticktext" => MONTH_NAMES, "showgrid" => false, "ticks" => ""),
        "yaxis" => axis("title" => titled("rain so far (mm)"), "rangemode" => "tozero")))
end

"""
Each year's value of variable `v`: for temperature and wind the year's mean (line) with the
band of its 5th to 95th percentile days and its lowest and highest days (dots); for rain the
year's total (bars) and its rain days. Shows how the climate changes over the years.
"""
function years_chart(df, v, meta, th)
    ys = sort(unique(df.year))
    by = [Float64.(df[df.year .== y, v]) for y ∈ ys]
    if v == "rain"
        tot = [sum(x) for x ∈ by]
        days = [count(>=(th.dry_day), x) for x ∈ by]
        return Dict("data" => [
                Dict("type" => "bar", "x" => ys, "y" => tot, "name" => "rainfall (mm)", "marker" => Dict("color" => rgba(meta.color, 0.75)),
                     "hovertemplate" => "%{x}: %{y:.0f} mm<extra></extra>"),
                Dict("x" => ys, "y" => days, "mode" => "lines+markers", "name" => "rain days (≥ $(th.dry_day) mm)", "yaxis" => "y2",
                     "line" => Dict("color" => C_DRY, "width" => 2), "marker" => Dict("color" => C_DRY, "size" => 6),
                     "hovertemplate" => "%{x}: %{y} rain days<extra></extra>")],
            "layout" => base_layout("height" => 320, "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.12),
                "margin" => Dict("l" => 60, "r" => 56, "t" => 26, "b" => 36), "bargap" => 0.25,
                "xaxis" => axis("tickformat" => "d", "dtick" => max(1, ceil(Int, length(ys) / 10)), "showgrid" => false),
                "yaxis" => axis("title" => titled("rainfall (mm)"), "rangemode" => "tozero"),
                "yaxis2" => axis("title" => titled("rain days"), "overlaying" => "y", "side" => "right", "showgrid" => false, "rangemode" => "tozero", "tickformat" => "d")))
    end
    q(p) = [quantile(x, p) for x ∈ by]
    u = meta.unit
    Dict("data" => [
            Dict("x" => ys, "y" => q(0.05), "mode" => "lines", "line" => Dict("width" => 0, "color" => meta.color), "hoverinfo" => "skip", "showlegend" => false),
            Dict("x" => ys, "y" => q(0.95), "mode" => "lines", "line" => Dict("width" => 0, "color" => meta.color), "fill" => "tonexty",
                 "fillcolor" => rgba(meta.color, 0.15), "name" => "5th–95th percentile", "hoverinfo" => "skip"),
            Dict("x" => ys, "y" => [mean(x) for x ∈ by], "mode" => "lines+markers", "name" => "mean",
                 "line" => Dict("color" => meta.color, "width" => 2.5), "marker" => Dict("color" => meta.color, "size" => 6),
                 "text" => ["5th–95th: $(fmtn(a)) to $(fmtn(b)) $u" for (a, b) ∈ zip(q(0.05), q(0.95))],
                 "hovertemplate" => "%{x}: mean %{y:.2f} $u; %{text}<extra></extra>"),
            Dict("x" => ys, "y" => [minimum(x) for x ∈ by], "mode" => "markers", "name" => "lowest / highest day",
                 "marker" => Dict("color" => meta.color, "size" => 5, "symbol" => "line-ew-open", "line" => Dict("width" => 2, "color" => meta.color)),
                 "hovertemplate" => "%{x}: lowest %{y:.1f} $u<extra></extra>"),
            Dict("x" => ys, "y" => [maximum(x) for x ∈ by], "mode" => "markers", "showlegend" => false,
                 "marker" => Dict("color" => meta.color, "size" => 5, "symbol" => "line-ew-open", "line" => Dict("width" => 2, "color" => meta.color)),
                 "hovertemplate" => "%{x}: highest %{y:.1f} $u<extra></extra>")],
         "layout" => base_layout("height" => 320, "legend" => Dict("orientation" => "h", "x" => 0, "y" => 1.12),
            "margin" => Dict("l" => 60, "r" => 12, "t" => 26, "b" => 36),
            "xaxis" => axis("tickformat" => "d", "dtick" => max(1, ceil(Int, length(ys) / 10)), "showgrid" => false),
            "yaxis" => axis("title" => titled("$(lowercase(meta.label)) ($u)"))))
end

# ---- the page ---------------------------------------------------------------------------------

"""
Everything the Explore tab shows for variable `v` of site `s`, for year `year`
(0: all years), with the thresholds `th` = (hot, dry_day, dry_month).
"""
function sim_page(s, v, year, th)
    df = sim_weather(s)
    isnothing(df) && return nothing
    years = sort(unique(df.year))
    d = year == 0 ? df : df[df.year .== year, :]
    nrow(d) == 0 && (d = df)
    meta = SW_VARS[v]
    has(c) = c in names(d)
    charts = Dict{String, Any}[]
    stats = Dict{String, Any}[]
    chances = Dict{String, Any}[]
    card(key, title, sub, chart) = push!(charts, Dict{String, Any}("key" => key, "title" => title, "sub" => sub, "chart" => chart))
    stat(l, x) = push!(stats, Dict{String, Any}("label" => l, "value" => x))
    if !has(v)
        return (years=years, charts=charts, stats=stats, chances=chances, extremes=extremes(df, years, th), note="The simulated weather has no $v.")
    end
    x = Float64.(d[!, v])
    if year == 0 && length(years) > 1
        card("years", "Over the years", v == "rain" ? "Each year's rainfall and rain days (rain ≥ $(th.dry_day) mm): how the climate changes from the first to the last year." :
             "Each year's mean, the band of its 5th to 95th percentile days, and its lowest and highest days: how the climate changes from the first to the last year.",
             years_chart(df, v, meta, th))
    end
    if v == "rain"
        wetd = x .>= th.dry_day
        amounts = x[wetd .& (x .> 0)]
        stat("wet days", "$(round(100 * mean(wetd); digits=1))% of days")
        if !isempty(amounts)
            stat("mean on a wet day", "$(fmtn(mean(amounts); d=1)) mm")
            stat("median on a wet day", "$(fmtn(quantile(amounts, 0.5); d=1)) mm")
            stat("95th percentile", "$(fmtn(quantile(amounts, 0.95); d=1)) mm")
            stat("largest", "$(fmtn(maximum(amounts); d=1)) mm")
            card("dist", "Rain on wet days", "Daily amounts of the wet days (rain ≥ $(th.dry_day) mm), on a log scale.",
                 dist_chart(amounts, meta.color, "rain on a wet day (mm, log scale)"; logx=true))
        end
        wa = [x[(d.month .== m) .& wetd .& (x .> 0)] for m ∈ 1:12]
        card("box", "Wet-day rain through the year", "Each month's wet-day amounts (rain ≥ $(th.dry_day) mm), on a log scale: the violin's width shows how often each amount occurs; the small box is the middle half, its line the median.",
             month_violin_chart(wa, meta.color, "rain on a wet day (mm, log scale)"; logy=true))
        cw, cd = wet_chances(d, wetd)
        card("chances", "Wet-day chances through the year",
             "Each month's pww and pwd in the generated days" * (year == 0 && length(years) > 1 ? " (line: mean over the years; band: their range)." : "."),
             chances_chart(cw, cd))
        mw, md = nanmean_row(cw), nanmean_row(cd)
        chances = map(1:12) do m
            w, p = mw[m], md[m]
            ok(v) = !isnan(v)
            Dict{String, Any}("month" => MONTH_NAMES[m], "pww" => ok(w) ? fmtn(w) : "–", "pwd" => ok(p) ? fmtn(p) : "–",
                "share" => ok(w) && ok(p) ? "$(round(Int, 100 * wet_fraction(w, p)))%" : "–",
                "wetspell" => ok(w) && w < 1 ? fmtn(wet_spell(w); d=1) : "–",
                "dryspell" => ok(p) && p > 0 ? fmtn(dry_spell(p); d=1) : "–")
        end
        ws, ds = runs(wetd), runs(.!wetd)
        stat("longest wet spell", "$(isempty(ws) ? 0 : maximum(ws)) days")
        stat("longest dry spell", "$(isempty(ds) ? 0 : maximum(ds)) days")
        card("spells", "Wet and dry spells", "How often runs of wet and of dry days of each length occur (wet day: rain ≥ $(th.dry_day) mm).",
             spell_chart(ws, ds))
        card("cum", "Rain through the year", length(years) < 5 ? "Rain accumulated from 1 January, one line per year." :
             "Rain accumulated from 1 January: each year a faint line, with the median of the years and the band of their 10th to 90th percentile" *
             (year == 0 ? "." : "; $year in orange."),
             cumrain_chart(df, length(years) < 5 && year != 0 ? [year] : years, year))
    else
        stat("mean", "$(fmtn(mean(x))) $(meta.unit)")
        stat("sd", "$(fmtn(std(x))) $(meta.unit)")
        stat("skew", fmtn(skew3(x)))
        stat("5th / 50th / 95th percentile", join(fmtn.(quantile(x, [0.05, 0.5, 0.95])), " / ") * " $(meta.unit)")
        stat("lowest / highest", "$(fmtn(minimum(x))) / $(fmtn(maximum(x))) $(meta.unit)")
        stat("day-to-day correlation (lag 1)", fmtn(lag1(x)))
        card("dist", "Distribution", "All days$(year == 0 ? "" : " of $year"): histogram and its smoothed curve.",
             dist_chart(x, meta.color, "$(lowercase(meta.label)) ($(meta.unit))"))
        card("box", "Through the year", "Each month's days: the violin's width shows how often each value occurs; the small box is the middle half, its line the median.",
             month_violin_chart([x[d.month .== m] for m ∈ 1:12], meta.color, meta.unit))
        card("lag", "Today against yesterday", "Each day against the day before: the tighter the cloud along the line, the stronger the day-to-day persistence (rlag).",
             lag_chart(x, meta.color, meta.unit))
        if v ∈ ("tmin", "tmax") && has("tmin") && has("tmax")
            rng = d.tmax .- d.tmin
            stat("daily range (tmax − tmin)", "mean $(fmtn(mean(rng))), smallest $(fmtn(minimum(rng))) °C")
        end
        if v == "tmax" && has("rain")
            wet = (d.rain .>= th.dry_day) .& (d.rain .> 0)
            if any(wet) && !all(wet)
                stat("mean tmax: wet / dry days", "$(fmtn(mean(d.tmax[wet]))) / $(fmtn(mean(d.tmax[.!wet]))) °C")
                card("train", "Monthly Tmax against monthly rain", "One dot per month of each year: the month's mean tmax against its total rainfall, with the least-squares line and the correlation r.",
                     tmax_rain_chart(d))
            end
        end
    end
    (years=years, charts=charts, stats=stats, chances=chances, extremes=extremes(df, years, th), note="")
end

"The extremes of each year: hot days, largest daily rain, longest dry spell, dry months, total rain."
function extremes(df, years, th)
    map(years) do y
        g = df[df.year .== y, :]
        r = Dict{String, Any}("year" => y)
        "tmax" in names(g) && (r["hot"] = count(>=(th.hot), g.tmax); r["tmax"] = fmtn(maximum(g.tmax); d=1))
        "tmin" in names(g) && (r["tmin"] = fmtn(minimum(g.tmin); d=1))
        "wind" in names(g) && (r["wind"] = fmtn(maximum(g.wind); d=1))
        if "rain" in names(g)
            ds = runs(g.rain .< th.dry_day)
            r["rainmax"] = fmtn(maximum(g.rain); d=1)
            r["dryspell"] = isempty(ds) ? 0 : maximum(ds)
            # only the months the year has (a partial year's missing months are not dry)
            r["drymonths"] = count(m -> any(g.month .== m) && sum(g.rain[g.month .== m]) < th.dry_month, 1:12)
            r["total"] = round(Int, sum(g.rain))
        end
        r
    end
end
