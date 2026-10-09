### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a2a81510-a24d-481d-a7de-a5bccbc77218
begin
    using Dates, Distributions, HypertextLiteral, PlutoUI
    using Random, SpecialFunctions, Statistics
end

# ╔═╡ 695aa061-f11d-4895-a66b-093f79ebdb2e
md"""
# Weather designer for MsiaGen

Design the daily weather statistics of a site (min. and max. air temperature, wind speed, and rainfall) for a **start year**, say how they **change by the end year**, then click **Generate** at the bottom to get the statistics file for MsiaGen.

Each part (min. and max. air temperature, wind speed, and rainfall) works in three layers:

1. **Set start year**: the whole year, a seasonal curve plus values that apply to all months.
2. **Each month**: change any month away from the whole-year values.
3. **Set last year**: the total change from the start year to the end year, spread evenly over the years between and applied to all months.

To modify an existing weather record instead of starting from typical Malaysian values, load a daily weather file under **Site Info**: every slider then starts from that file's average year.

The points and rain days in the charts are **illustrative**: drawn with MsiaGen's own equations so you can see what the numbers mean. They are not the generated weather. Use the table of contents on the right to jump between parts.
"""

# ╔═╡ 20609a8a-99cc-4411-b1fa-8f1e2a709609
TableOfContents()

# ╔═╡ c6cfd927-288b-4ada-b19c-5dc6fccc2a1c
md"""
**Start from a daily weather file (optional):** $(@bind wfile FilePicker([MIME("text/csv"), MIME("text/plain")]))

Any daily weather CSV with columns `year`, `month`, `day` and some of `tmin`, `tmax`, `wind`, `rain`, such as a site's `<site>-obs.csv` or `<site>-sim.csv`. Its years are averaged into one typical year, and every slider below starts from it. Site, latitude and years above are filled in from the file.
"""

# ╔═╡ a4424332-1266-4ad7-8559-e0ba6ceb7a54
begin
    const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    # change from the start value, shown with its sign: +0.1, −0.1, 0
    fmt_change(v) = v == 0 ? "0" : (v > 0 ? "+" : "−") * string(abs(v))
    # slider steps from lo to hi (adding 0.0 turns a rounded -0.0 into 0.0)
    steps(lo, hi, step) = round.(range(lo, hi; step=step); digits=3) .+ 0.0
    # days in each month of year y, and the month of every day
    function calendar(y)
        mdays = [daysinmonth(y, m) for m ∈ 1:12]
        (mdays=mdays, ndays=sum(mdays),
         month_of=reduce(vcat, [fill(m, mdays[m]) for m ∈ 1:12]))
    end
    # chart colours (reference palette of the dataviz guidelines)
    const BLUE = "#2a78d6"
    const DARK_BLUE = "#104281"
    const INK = "#0b0b0b"
    const INK2 = "#52514e"
    const GRID = "#e1e0d9"
    const MUTED = "#898781"
    const SURFACE = "#fcfcfb"
    const RED = "#e34948"
    const DARK_RED = "#a32525"
    # each temperature keeps its colour in every chart: (points and band, mean)
    const T_COLORS = Dict("tmin" => (BLUE, DARK_BLUE), "tmax" => (RED, DARK_RED))
    # How a value moves from its start-year value to its end-year target:
    # g(f) is the share of the change reached at fraction f of the way from
    # the start year (f = 0) to the end year (f = 1)
    const PATHS = ["linear" => "Linear", "accel" => "Accelerating",
                   "level" => "Levelling off", "scurve" => "S-curve"]
    path_g(path, f) = path == "accel" ? f^2 :
                      path == "level" ? 1 - (1 - f)^2 :
                      path == "scurve" ? 3f^2 - 2f^3 : f
    toward(a, b, g) = a + g * (b - a)
    # a day of the year moves the shorter way round the year (350 -> 10 is +25)
    toward_day(a, b, g) = a + g * (mod(b - a + 182, 365) - 182)
    # a day of the year as shown: within 1-365
    show_day(d) = mod1(round(Int, d), 365)

    # results tables: centre every column
    const CENTRED = @htl("""<style>
        table.centred th, table.centred td { text-align: center !important; }
        </style>""")

    # Charts are SVG, drawn by the browser, so no plotting package is needed
    # (a plotting package took most of the notebook's start-up time).
    svg_text(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
    r1(x) = round(x; digits=1)
    # tidy axis ticks from lo to hi
    function nice_ticks(lo, hi; n=5)
        raw = (hi - lo) / n
        mag = 10.0^floor(log10(raw))
        f = raw / mag
        step = mag * (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10)
        [round(k * step; digits=6) for k ∈ ceil(Int, lo / step - 1e-9):floor(Int, hi / step + 1e-9)]
    end
    tick_label(v) = isinteger(v) ? string(Int(v)) : string(v)
    # an SVG chart of w × h units, scaled to the cell's width
    svg(body, w, h) = HTML("""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $w $h"
        width="100%" style="max-width: $(w)px; height: auto; display: block;
        background: $SURFACE; font-family: system-ui, -apple-system, 'Segoe UI',
        Helvetica, sans-serif;">$body</svg>""")
    # A legend starting at (x, y): entries (kind, colour, label), with kind
    # :dot, :line, :dash, :dotted, :band (light fill) or :box; entries that
    # would pass `maxx` go on a new row
    function svg_legend(io, x, y, entries; title=nothing, maxx=Inf)
        x0 = x
        if !isnothing(title)
            print(io, """<text x="$x" y="$(y + 5)" font-size="14" font-weight="600" fill="$INK2">$(svg_text(title))</text>""")
            x += 9 * length(title) + 20
        end
        for (kind, col, label) ∈ entries
            w = 26 + 7.5 * length(label)
            if x > x0 && x + w > maxx
                x, y = x0, y + 24
            end
            if kind == :dot
                print(io, """<circle cx="$(x + 8)" cy="$y" r="4" fill="$col"/>""")
            elseif kind ∈ (:box, :band)
                fill = kind == :band ? """fill="$col" fill-opacity="0.12\"""" :
                       """fill="$col" stroke="$GRID" stroke-width="0.5\""""
                print(io, """<rect x="$x" y="$(y - 7)" width="16" height="14" $fill/>""")
            else
                dash = kind == :dash ? "6 4" : kind == :dotted ? "2 3" : "none"
                print(io, """<line x1="$x" x2="$(x + 20)" y1="$y" y2="$y" stroke="$col" stroke-width="2.5" stroke-dasharray="$dash"/>""")
            end
            print(io, """<text x="$(x + 26)" y="$(y + 5)" font-size="14" fill="$INK2">$(svg_text(label))</text>""")
            x += w + 20
        end
    end

    # Chart of one year's illustrative days `pts` around its monthly
    # statistics `s` (μ, σ, mdays, ndays): the ±1 sd band and mean of each
    # month; optionally `ref`, another year's monthly means (dashed),
    # `floor_at`, a lower limit (dotted; bands stop there), and `other`, the
    # monthly statistics of another variable (`other.s`, named `other.label`,
    # in `other.color`), drawn faded behind as its mean and ±1 sd band.
    # `colors`: (points and band, mean) of the shown variable; `label`: its
    # name in the legend. Hover a day or the other variable's band to see
    # its value.
    function day_chart(s, pts, title, ylabel, unit; ref=nothing, floor_at=nothing,
                       floor_label="", other=nothing, colors=(BLUE, DARK_BLUE), label="")
        col, dark = colors
        W, H = 1000, isnothing(other) ? 470 : 590
        L, R, T, B = 70, 20, 44, 88            # margins
        mdays, ndays = s.mdays, s.ndays
        starts = cumsum([1; mdays[1:end-1]])
        band_lo = isnothing(floor_at) ? s.μ .- s.σ : max.(floor_at, s.μ .- s.σ)
        lo, hi = extrema([pts; band_lo; s.μ .+ s.σ; something(ref, Float64[]);
                          isnothing(other) ? Float64[] : [other.s.μ .- other.s.σ; other.s.μ .+ other.s.σ]])
        pad = 0.05 * max(hi - lo, 1e-6)
        lo, hi = isnothing(floor_at) ? lo - pad : 0.0, hi + pad
        X(d) = L + (d - 0.5) / ndays * (W - L - R)
        Y(v) = T + (hi - v) / (hi - lo) * (H - T - B)
        io = IOBuffer()
        print(io, """<text x="$L" y="26" font-size="17" font-weight="600" fill="$INK">$(svg_text(title))</text>""")
        for v ∈ nice_ticks(lo, hi)
            y = r1(Y(v))
            print(io, """<line x1="$L" x2="$(W - R)" y1="$y" y2="$y" stroke="$GRID"/>""",
                      """<text x="$(L - 8)" y="$(y + 5)" font-size="14" fill="$INK2" text-anchor="end">$(tick_label(v))</text>""")
        end
        print(io, """<text transform="translate(20 $(r1((T + H - B) / 2))) rotate(-90)" font-size="14" fill="$INK2" text-anchor="middle">$(svg_text(ylabel))</text>""",
                  """<line x1="$L" x2="$L" y1="$T" y2="$(H - B)" stroke="$GRID"/>""",
                  """<line x1="$L" x2="$(W - R)" y1="$(H - B)" y2="$(H - B)" stroke="$GRID"/>""")
        for m ∈ 1:12
            x0, x1 = r1(X(starts[m] - 0.5)), r1(X(starts[m] + mdays[m] - 0.5))
            m > 1 && print(io, """<line x1="$x0" x2="$x0" y1="$T" y2="$(H - B)" stroke="$GRID"/>""")
            print(io, """<text x="$(r1((x0 + x1) / 2))" y="$(H - B + 22)" font-size="14" fill="$INK2" text-anchor="middle">$(MONTHS[m])</text>""")
            if !isnothing(other)
                o = other.s
                ot, ob, om = r1(Y(o.μ[m] + o.σ[m])), r1(Y(o.μ[m] - o.σ[m])), r1(Y(o.μ[m]))
                tip = "$(other.label), $(MONTHS[m]): mean $(round(o.μ[m]; digits=2)) $unit, sd $(o.σ[m]) $unit"
                print(io, """<rect x="$x0" y="$ot" width="$(r1(x1 - x0))" height="$(r1(ob - ot))" fill="$(other.color)" fill-opacity="0.12"><title>$(svg_text(tip))</title></rect>""",
                          """<line x1="$x0" x2="$x1" y1="$om" y2="$om" stroke="$(other.color)" stroke-opacity="0.7" stroke-width="2"/>""")
            end
            yt, yb, ym = r1(Y(s.μ[m] + s.σ[m])), r1(Y(band_lo[m])), r1(Y(s.μ[m]))
            print(io, """<rect x="$x0" y="$yt" width="$(r1(x1 - x0))" height="$(r1(yb - yt))" fill="$col" fill-opacity="0.12"/>""",
                      """<line x1="$x0" x2="$x1" y1="$ym" y2="$ym" stroke="$dark" stroke-width="2.5"/>""")
            if !isnothing(ref)
                yr = r1(Y(ref[m]))
                print(io, """<line x1="$x0" x2="$x1" y1="$yr" y2="$yr" stroke="$MUTED" stroke-width="2" stroke-dasharray="6 4"/>""")
            end
        end
        if !isnothing(floor_at)
            yf = r1(Y(floor_at))
            print(io, """<line x1="$L" x2="$(W - R)" y1="$yf" y2="$yf" stroke="$MUTED" stroke-dasharray="2 3"/>""")
        end
        coords = join(("$(r1(X(d))),$(r1(Y(v)))" for (d, v) ∈ enumerate(pts)), " ")
        print(io, """<polyline points="$coords" fill="none" stroke="$col" stroke-opacity="0.35" stroke-width="0.8"/>""")
        for m ∈ 1:12, k ∈ 1:mdays[m]
            d = starts[m] + k - 1
            print(io, """<circle cx="$(r1(X(d)))" cy="$(r1(Y(pts[d])))" r="2.5" fill="$col"><title>$(MONTHS[m]) $k: $(round(pts[d]; digits=2)) $unit</title></circle>""")
        end
        pre = isempty(label) ? "" : "$(label) "
        entries = [(:dot, col, "illustrative $(pre)day"), (:line, dark, "$(pre)monthly mean"),
                   (:band, col, "$(pre)±1 sd")]
        isnothing(floor_at) || push!(entries, (:dotted, MUTED, floor_label))
        isnothing(ref) || push!(entries, (:dash, MUTED, "start-year $(pre)monthly mean"))
        svg_legend(io, L, H - B + 50, entries; maxx=W - R)
        # the other variable on a row of its own
        isnothing(other) || svg_legend(io, L, H - B + 74, [(:line, other.color, "$(other.label) monthly mean"),
                                                         (:band, other.color, "$(other.label) ±1 sd")])
        svg(String(take!(io)), W, H)
    end
end;

# ╔═╡ 1f5fb93a-a616-4cd2-8eb4-2d68e9f6a07e
begin
    # Reads a daily weather CSV: an optional latitude line, '#' comment lines,
    # a header, then numbers (missing values become NaN).
    function parse_weather(text::AbstractString)
        lines = [strip(l) for l ∈ split(text, '\n')]
        lines = [l for l ∈ lines if !isempty(l) && !startswith(l, "#")]
        lat = tryparse(Float64, first(lines))
        isnothing(lat) || (lines = lines[2:end])
        hdr = lowercase.(strip.(split(first(lines), ',')))
        vals = [[something(tryparse(Float64, strip(v)), NaN) for v ∈ split(l, ',')]
                for l ∈ lines[2:end]]
        cols = Dict(h => [length(r) >= i ? r[i] : NaN for r ∈ vals] for (i, h) ∈ enumerate(hdr))
        all(haskey(cols, k) for k ∈ ("year", "month", "day")) ||
            error("the file needs columns year, month and day")
        (lat=lat, cols=cols)
    end

    # Monthly statistics, computed as MsiaGen computes them (input.jl, utils.jl)
    nanmean(x) = (y = filter(!isnan, x); isempty(y) ? NaN : mean(y))
    function lag1(x)
        m = mean(x)
        d = sum((x .- m) .^ 2)
        d ≈ 0 ? 0.0 : sum((x[1:end-1] .- m) .* (x[2:end] .- m)) / d
    end
    function skew3(x)
        m = mean(x)
        m2 = mean((x .- m) .^ 2)
        m2 ≈ 0 ? 0.0 : mean((x .- m) .^ 3) / m2^1.5
    end
    # pww and pwd over the days that have a next day (all but the last);
    # no dry days (every day wet): pwd = 1, as in MsiaGen
    function rain_probs(r)
        wet = r .> 0
        prev, next = wet[1:end-1], wet[2:end]
        nw, nd = count(prev), count(.!prev)
        (nw > 0 ? count(prev .& next) / nw : 0.0, nd > 0 ? count(.!prev .& next) / nd : 1.0)
    end

    # The file's average year: each statistic of each month, from every year,
    # averaged over the years. A variable missing from the file is `nothing`.
    function average_year(f)
        c = f.cols
        order = sortperm(collect(zip(c["year"], c["month"], c["day"])))
        yrs = sort(unique(filter(!isnan, c["year"])))
        isempty(yrs) && error("no years found in the file")
        month_vals(col, y, m) = filter(!isnan, c[col][[i for i ∈ order
                                       if c["year"][i] == y && c["month"][i] == m]])
        per_month(col, stat) =
            [nanmean([(x = month_vals(col, y, m); length(x) >= 3 ? stat(x) : NaN) for y ∈ yrs])
             for m ∈ 1:12]
        sd(x) = max(0.01, std(x))
        temp(col) = haskey(c, col) ?
            (μ=per_month(col, mean), σ=per_month(col, sd), ρ=per_month(col, lag1),
             γ=per_month(col, skew3)) : nothing
        # a variable with a month that has no data in any year is left out
        complete(x) = isnothing(x) || all(v -> !any(isnan, v), values(x)) ? x : nothing
        (lat=f.lat, years=Int.(yrs), tmin=complete(temp("tmin")), tmax=complete(temp("tmax")),
         wind=complete(haskey(c, "wind") ?
             (μ=per_month("wind", mean), σ=per_month("wind", sd), ρ=per_month("wind", lag1)) : nothing),
         rain=complete(haskey(c, "rain") ?
             (total=per_month("rain", sum), pww=per_month("rain", x -> rain_probs(x)[1]),
              pwd=per_month("rain", x -> rain_probs(x)[2])) : nothing))
    end

    # Fitting the designer's seasonal curves to 12 monthly values, on a
    # 365-day year (so the start positions do not depend on the start year)
    snap(r, x) = r[argmin(abs.(r .- x))]
    const CAL365 = calendar(2001)
    const MONTH_COS = [mean(cos.(2π .* findall(CAL365.month_of .== m) ./ 365)) for m ∈ 1:12]
    const MONTH_SIN = [mean(sin.(2π .* findall(CAL365.month_of .== m) ./ 365)) for m ∈ 1:12]
    day_angle(day) = 2π * day / 365
    angle_day(θ) = mod(round(Int, θ * 365 / (2π)) - 1, 365) + 1

    # monthly means of the curve M + A cos(2π (day - peak) / 365), and the
    # curve that fits 12 monthly means best
    curve_months(M, A, peak) = M .+ A .* (cos(day_angle(peak)) .* MONTH_COS .+
                                          sin(day_angle(peak)) .* MONTH_SIN)
    function fit_curve(μ)
        β = [ones(12) MONTH_COS MONTH_SIN] \ μ
        (M=β[1], A=hypot(β[2], β[3]), peak=angle_day(atan(β[3], β[2])))
    end
    # monthly shares of the rain curve 1 - s cos(2π (day - dry) / 365), and
    # the curve that fits 12 monthly totals best
    function rain_shares(seas, dry)
        curve = 1 .- (seas / 100) .* cos.(2π .* ((1:365) .- dry) ./ 365)
        [sum(curve[CAL365.month_of .== m]) for m ∈ 1:12] ./ sum(curve)
    end
    function fit_rain(total)
        β = [ones(12) MONTH_COS MONTH_SIN] \ (total ./ CAL365.mdays)
        a, b = β[2] / β[1], β[3] / β[1]
        (seas=100 * hypot(a, b), dry=angle_day(atan(-b, -a)))
    end
end;

# ╔═╡ d49f9ae6-9b07-484d-a105-13b3fbac522f
# The loaded file's average year (or nothing, or the error when it cannot
# be read)
baseline = isnothing(wfile) ? nothing :
           try average_year(parse_weather(String(copy(wfile["data"])))) catch err err end;

# ╔═╡ 6ad69763-051b-4d69-a99f-69c495577cab
# Defaults of the Site Info fields: from the loaded weather file (its name
# without "-obs", "-sim" or "-stats", its latitude line, and its first and
# last years), otherwise empty site, 3.0° and 2000-2020
site_defaults = let b = baseline isa NamedTuple ? baseline : nothing
    if isnothing(b)
        (site="", lat=3.0, start=2000, stop=2020)
    else
        name = replace(replace(wfile["name"], r"\.[^.]*$" => ""), r"-(obs|sim|stats)$"i => "")
        (site=name,
         lat=isnothing(b.lat) ? 3.0 : clamp(round(b.lat; digits=3), -10.0, 10.0),
         start=clamp(first(b.years), 1950, 2100), stop=clamp(last(b.years), 1950, 2100))
    end
end;

# ╔═╡ 86792682-a96d-4ab2-a321-4199e2e5c0b1
md"""
# Site Info

**Site:** $(@bind site TextField(; default=site_defaults.site, placeholder="e.g. Kluang"))
  **Latitude (°):** $(@bind lat NumberField(-10.0:0.001:10.0; default=site_defaults.lat))

  **Start year:** $(@bind start_year NumberField(1950:2100; default=site_defaults.start))
  **End year:** $(@bind end_year NumberField(1950:2100; default=site_defaults.stop))
"""

# ╔═╡ 1db20637-6907-44cc-8728-400469155d76
if isnothing(baseline)
    md"_No file loaded: sliders start from typical Malaysian values._"
elseif baseline isa Exception
    Markdown.parse("""
    !!! danger "Could not read $(wfile["name"])"
        $(sprint(showerror, baseline))
    """)
else
    let vars = ("tmin", "tmax", "wind", "rain")
        found = [v for v ∈ vars if !isnothing(getproperty(baseline, Symbol(v)))]
        absent = setdiff(collect(vars), found)
        note_absent = isempty(absent) ? "" : " ($(join(absent, ", ")) missing or incomplete in the file: typical values used)"
        note_lat = " Site, latitude and years above were filled in from the file" *
                   (isnothing(baseline.lat) ? " (it has no latitude line, so check the latitude)." : ".")
        Markdown.parse("""
        !!! info "Loaded $(wfile["name"])"
            Average of $(length(baseline.years)) year(s), $(first(baseline.years))–$(last(baseline.years)). Sliders start from it for $(join(found, ", "))$(note_absent).$(note_lat)
        """)
    end
end

# ╔═╡ 2ac2e9cc-c555-46bc-9bf0-47bc6340192d
# years to generate, and how far year y is from the start year to the end
# year (0 in the start year, 1 in the end year)
begin
    years = start_year:max(start_year, end_year)
    frac(y) = length(years) == 1 ? 0.0 : (y - start_year) / (last(years) - start_year)
end;

# ╔═╡ cf725ce2-bebe-4b48-8405-d369e1cef355
md"""
# 1. Air temperature

MsiaGen needs four statistics per month for each of minimum (`tmin`) and maximum (`tmax`) air temperature: **mean**, **sd**, **rlag** and **skew**. Choose **Tmax** or **Tmin** below to see and set its sliders; the settings of both are kept.

Start values when no weather file is loaded (`sd`, `rlag` and `skew` are typical of Malaysian sites):

| | annual mean (°C) | amplitude (°C) | warmest day | sd (°C) | rlag | skew |
|---|---|---|---|---|---|---|
| tmin | 22 | 1 | 183 (mid-year) | 0.7 | 0.25 | −0.1 |
| tmax | 33 | 1 | 183 (mid-year) | 1.25 | 0.2 | −0.8 |

- `sd`: spread of days around the monthly mean (the shaded band is ±1 sd)
- `rlag`: how much a day follows the day before (high = smooth trails, low = jumpy)
- `skew`: lopsidedness (negative = more cool dips, positive = more hot spikes)
"""

# ╔═╡ 7a6dc0d6-5d28-4ea3-a728-3a873958d341
@htl("""<p><b>Select air temperature type:</b> <span id="tvar-select">$(@bind tvar Select(["tmin" => "tmin", "tmax" => "tmax"]))</span></p>""")

# ╔═╡ a28b530c-ac28-4989-92aa-448efae7e9ad
begin
    const TVARS = ["tmin", "tmax"]      # same order as the dropdown
    # Start values of the whole year, used when no weather file is loaded.
    # sd, rlag and skew are the medians of all monthly values in the 23
    # sites' stats files.
    const T_START = Dict("tmax" => (M=33.0, A=1.0, peak=183, sd=1.25, rlag=0.2, skew=-0.8),
                         "tmin" => (M=22.0, A=1.0, peak=183, sd=0.7, rlag=0.25, skew=-0.1))
    # valid limits of each statistic in MsiaGen
    const T_LIMITS = (sd=(0.1, 4.0), rlag=(0.0, 0.95), skew=(-2.0, 2.0))
    # slider values: whole year, and month changes (d...)
    const T_RANGE = (M=15.0:0.1:40.0, A=0.0:0.05:5.0, peak=1:365,
                     sd=steps(T_LIMITS.sd..., 0.05), rlag=steps(T_LIMITS.rlag..., 0.01),
                     skew=steps(T_LIMITS.skew..., 0.05),
                     dm=steps(-4.0, 4.0, 0.1), dsd=steps(-2.0, 2.0, 0.05),
                     drl=steps(-1, 1, 0.1), dsk=steps(-2.0, 2.0, 0.05))
    # Shows only the sliders of the variable chosen in the dropdown (tmax or
    # tmin). Both sets stay on the page, so their settings are kept.
    const T_TOGGLE = @htl("""<script>
        const root = currentScript.parentElement;
        const find = () => document.querySelector("#tvar-select select");
        function update() {
            const sel = find();
            if (!sel) return;
            root.querySelectorAll(".tvar-part").forEach(d =>
                d.style.display = (d.dataset.var == sel.selectedIndex) ? "" : "none");
        }
        function attach() {
            const sel = find();
            if (!sel) { setTimeout(attach, 200); return; }
            sel.addEventListener("input", update);
            invalidation.then(() => sel.removeEventListener("input", update));
            update();
        }
        attach();
        </script>""")
end

# ╔═╡ 6d7815aa-155b-438f-999d-f1b78b356d3a
# Start positions of the temperature sliders: from the loaded file's average
# year (whole-year curve fitted to its monthly means; each month slider holds
# that month's remainder, all snapped to the slider steps), or the start
# values with no month changes.
T_INIT = Dict(map(TVARS) do v
    b = baseline isa NamedTuple ? getproperty(baseline, Symbol(v)) : nothing
    if isnothing(b)
        d = T_START[v]
        v => (whole=d, dm=zeros(12), dsd=zeros(12), drl=zeros(12), dsk=zeros(12))
    else
        r = T_RANGE
        f = fit_curve(b.μ)
        M, A = snap(r.M, f.M), snap(r.A, f.A)
        sd, rl, sk = snap(r.sd, mean(b.σ)), snap(r.rlag, mean(b.ρ)), snap(r.skew, mean(b.γ))
        v => (whole=(M=M, A=A, peak=f.peak, sd=sd, rlag=rl, skew=sk),
              dm=snap.(Ref(r.dm), b.μ .- curve_months(M, A, f.peak)),
              dsd=snap.(Ref(r.dsd), b.σ .- sd), drl=snap.(Ref(r.drl), b.ρ .- rl),
              dsk=snap.(Ref(r.dsk), b.γ .- sk))
    end
end);

# ╔═╡ 4d2a0300-e697-4deb-8580-11e978893162
Markdown.parse("## a. ($(tvar)) Set start year ($(start_year))")

# ╔═╡ 05e5deb2-0c98-415d-a129-b8db81c24db4
@bind twhole PlutoUI.combine() do Child
    parts = map(enumerate(TVARS)) do (i, v)
        d, r = T_INIT[v].whole, T_RANGE
        c(k, widget) = Child(Symbol("$(v)_$(k)"), widget)
        @htl("""
        <div class="tvar-part" data-var="$(i - 1)">
        <table>
        <tr><td>Annual mean (°C)</td>
            <td>$(c(:M, PlutoUI.Slider(r.M; default=d.M, show_value=true)))</td></tr>
        <tr><td>Amplitude (°C)</td>
            <td>$(c(:A, PlutoUI.Slider(r.A; default=d.A, show_value=true)))</td></tr>
        <tr><td>Warmest day (day of year)</td>
            <td>$(c(:peak, PlutoUI.Slider(r.peak; default=d.peak, show_value=true)))</td></tr>
        <tr><td>sd (°C)</td>
            <td>$(c(:sd, PlutoUI.Slider(r.sd; default=d.sd, show_value=true)))</td></tr>
        <tr><td>rlag</td>
            <td>$(c(:rlag, PlutoUI.Slider(r.rlag; default=d.rlag, show_value=true)))</td></tr>
        <tr><td>skew</td>
            <td>$(c(:skew, PlutoUI.Slider(r.skew; default=d.skew, show_value=true)))</td></tr>
        </table>
        </div>
        """)
    end
    @htl("""<div>$(parts)$(T_TOGGLE)</div>""")
end

# ╔═╡ 192267b9-fadc-4bc5-9bb0-afa0712d1e49
@bind tnew CounterButton("New sample of points")

# ╔═╡ c0544ade-de13-48eb-937d-ef45ce85700c
md"""
## b. ($(tvar)) Each month
"""

# ╔═╡ 77485903-cd54-41ae-ad0f-460ecb1e6058
# Each slider is a change from the whole-year value of the start year.
# Results beyond T_LIMITS are capped (and reported in the output).
@bind tmonthly PlutoUI.combine() do Child
    parts = map(enumerate(TVARS)) do (i, v)
        d, r = T_INIT[v], T_RANGE
        c(k, m, widget) = Child(Symbol("$(v)_$(k)$(m)"), widget)
        rows = map(1:12) do m
            @htl("""
            <tr>
              <td><b>$(MONTHS[m])</b></td>
              <td>$(c("dm", m, PlutoUI.Slider(r.dm; default=d.dm[m], show_value=true)))</td>
              <td>$(c("sd", m, PlutoUI.Slider(r.dsd; default=d.dsd[m], show_value=true)))</td>
              <td>$(c("rl", m, PlutoUI.Slider(r.drl; default=d.drl[m], show_value=true)))</td>
              <td>$(c("sk", m, PlutoUI.Slider(r.dsk; default=d.dsk[m], show_value=true)))</td>
            </tr>
            """)
        end
        @htl("""
        <div class="tvar-part" data-var="$(i - 1)">
        <table>
        <tr><th>Month</th><th>mean change (°C)</th><th>sd change (°C)</th><th>rlag change</th><th>skew change</th></tr>
        $(rows)
        </table>
        </div>
        """)
    end
    @htl("""<div>$(parts)$(T_TOGGLE)</div>""")
end

# ╔═╡ 3cc149d9-614e-4ad1-a2e8-fe57e395a2aa
let e = max(start_year, end_year)
    Markdown.parse("""
    ## c. ($(tvar)) Set last year ($(e))

    Set the targets for $(e), and the path the values take from $(start_year) to get there; the years between follow that path, in all months. Targets start at the $(start_year) values, so changing section 1 resets them. The preview below shows $(e).
    """)
end

# ╔═╡ 4b6ac062-0db7-4d6d-8649-973e913c59d8
# fixed random numbers, so points move smoothly when a slider changes
t_draws = let rng = MersenneTwister(1234 + tnew)
    (u1=randn(rng, 366), u2=randn(rng, 366), uf=rand(rng, 366))
end;

# ╔═╡ eea66297-900f-48a8-af1d-4881cb5fd1cd
begin
    # Monthly temperature statistics of variable v (tmax or tmin) in year y.
    # The whole-year values move from the start year (section 1) towards the
    # end-year targets (section 3) along the chosen path; the year's seasonal
    # curve gives its 12 monthly means; the month sliders are then added.
    # Values beyond T_LIMITS are capped, and the year is flagged.
    # `whole`, `target` and `monthly` are the slider values of sections 1, 3
    # and 2.
    function t_year_stats(v, y, whole, target, monthly)
        w(k) = whole[Symbol("$(v)_$(k)")]
        t(k) = target[Symbol("$(v)_$(k)")]
        mc(k) = [monthly[Symbol("$(v)_$(k)$(m)")] for m ∈ 1:12]
        g = path_g(t(:path), frac(y))
        M = toward(w(:M), t(:M), g)
        A = max(0.0, toward(w(:A), t(:A), g))
        peak = toward_day(w(:peak), t(:peak), g)
        cal = calendar(y)
        curve = M .+ A .* cos.(2π .* ((1:cal.ndays) .- peak) ./ cal.ndays)
        μ = [mean(curve[cal.month_of .== m]) for m ∈ 1:12] .+ mc("dm")
        raw = (sd=toward(w(:sd), t(:sd), g) .+ mc("sd"),
               rlag=toward(w(:rlag), t(:rlag), g) .+ mc("rl"),
               skew=toward(w(:skew), t(:skew), g) .+ mc("sk"))
        capped(k) = round.(clamp.(raw[k], T_LIMITS[k]...); digits=3)
        (year=y, M=M, A=A, peak=peak, μ=μ,
         σ=capped(:sd), ρ=capped(:rlag), γ=capped(:skew),
         mdays=cal.mdays, ndays=cal.ndays,
         clamped=any(any(raw[k] .!= clamp.(raw[k], T_LIMITS[k]...)) for k ∈ keys(raw)))
    end
    
    # Skewed residuals with mean 0 and SD `sde`, using MsiaGen's two methods
    # (gentemp.jl: normal_skew and high_skew), driven by fixed random numbers.
    function skewed_residuals(sde, skew, u1, u2, uf)
        if abs(skew) < 0.995272          # skewed normal
            sk2_3 = abs(skew)^(2 / 3)
            n1 = 0.5 * π * sk2_3
            n2 = sk2_3 + ((4 - π) / 2)^(2 / 3)
            δ = copysign(sqrt(n1 / n2), skew)
            shape = δ / sqrt(1 - δ^2)
            dlt = shape / sqrt(1 + shape^2)
            scale = sqrt(sde^2 / (1 - 2 * dlt^2 / π))
            loc = -scale * sqrt(2 / π) * dlt
            x = copy(u1)
            flip = u2 .> shape .* u1
            x[flip] .*= -1
            loc .+ scale .* x
        else                             # F-distribution for high skews
            sk = abs(skew)
            df2 = 500
            d6, d4, d2 = df2 - 6, df2 - 4, df2 - 2
            a = sqrt(-32 * d4 + sk^2 * d6^2)
            df1 = -(d2 * (-d6 * sk + a)) / (2 * a)
            f = quantile.(FDist(df1, df2), uf)
            r = (f .- mean(f)) .* sde ./ std(f)
            skew < 0 && (r = -r)
            r .- mean(r)
        end
    end

    # One illustrative year of temperature: MsiaGen's lag-1 autoregression,
    # month by month (gentemp.jl: temp_dist!), with fixed random numbers.
    function t_sample_year(s, draws)
        x = zeros(s.ndays)
        yr_avg = sum(s.μ .* s.mdays) / s.ndays
        t0 = 1
        for m ∈ 1:12
            idx = t0:t0+s.mdays[m]-1
            sde = s.σ[m] * sqrt(1 - s.ρ[m]^2)
            c = s.μ[m] * (1 - s.ρ[m])
            e = skewed_residuals(sde, s.γ[m], draws.u1[idx], draws.u2[idx], draws.uf[idx])
            for (k, t) ∈ enumerate(idx)
                prev = (t == 1) ? yr_avg : x[t-1]
                x[t] = c + s.ρ[m] * prev + e[k]
            end
            t0 += s.mdays[m]
        end
        x
    end    

    # Chart of one year's illustrative days around its monthly statistics;
    # `ref`: monthly means of another year, drawn as a dashed reference
    # the other temperature variable: tmax for tmin, and tmin for tmax
    other_tvar(v) = v == "tmin" ? "tmax" : "tmin"
    # `other`: the other variable's statistics of the same year, drawn faded
    # behind. Tmin is always blue and Tmax red.
    t_chart(s, pts, v, title; ref=nothing, other=nothing) =
        day_chart(s, pts, title, "air temperature (°C)", "°C"; ref=ref,
                  colors=T_COLORS[v], label=v,
                  other=isnothing(other) ? nothing :
                        (s=other, label=other_tvar(v), color=first(T_COLORS[other_tvar(v)])))

    # Statistics of every year, for tmin and tmax
    all_t_stats(years, whole, target, monthly) =
        Dict(v => [t_year_stats(v, y, whole, target, monthly) for y ∈ years] for v ∈ TVARS)

    # Sliders of the end-year targets of the whole-year values and the path
    # to them (section 3), starting at the start-year values `whole`
    function t_target_widget(whole)
            PlutoUI.combine() do Child
            parts = map(enumerate(TVARS)) do (i, v)
                r = T_RANGE
                now(k) = whole[Symbol("$(v)_$(k)")]
                c(k, widget) = Child(Symbol("$(v)_$(k)"), widget)
                sl(k) = c(k, PlutoUI.Slider(r[k]; default=now(k), show_value=true))
                row(label, k) = @htl("""<tr><td>$(label)</td><td>$(sl(k))</td><td>$(now(k))</td></tr>""")
                @htl("""
                <div class="tvar-part" data-var="$(i - 1)">
                <p><b>Path:</b> $(c(:path, Select(PATHS)))</p>
                <table>
                <tr><th></th><th>Target in end year</th><th>Start year</th></tr>
                $(row("Annual mean (°C)", :M))
                $(row("Amplitude (°C)", :A))
                $(row("Warmest day (day of year)", :peak))
                </table>
                <details><summary>More settings: sd, rlag, skew</summary>
                <table>
                <tr><th></th><th>Target in end year</th><th>Start year</th></tr>
                $(row("sd (°C)", :sd))
                $(row("rlag", :rlag))
                $(row("skew", :skew))
                </table>
                </details>
                </div>
                """)
            end
            @htl("""<div>$(parts)$(T_TOGGLE)</div>""")
        end
    end

    # Skew-normal distribution with a month's mean, sd and skew, as MsiaGen
    # draws daily temperatures (|skew| taken as at most 0.99)
    function skewnormal(avg, sd, skew)
        sk = clamp(skew, -0.99, 0.99)
        sk2_3 = abs(sk)^(2 / 3)
        δ = copysign(sqrt(0.5π * sk2_3 / (sk2_3 + ((4 - π) / 2)^(2 / 3))), sk)
        ω = sd / sqrt(1 - 2δ^2 / π)
        SkewNormal(avg - ω * sqrt(2 / π) * δ, ω, δ / sqrt(1 - δ^2))
    end

    # Chance that a day's tmin reaches its tmax, from the month's (mean, sd,
    # skew) of tmax `tx` and tmin `tn`, drawn independently as MsiaGen does.
    # Within about 15% of the days MsiaGen repairs, at 5 observed sites.
    function p_cross(tx, tn)
        dx, dn = skewnormal(tx...), skewnormal(tn...)
        w = 8 * max(tx[2], tn[2])
        t = range(min(tx[1], tn[1]) - w, max(tx[1], tn[1]) + w; length=2001)
        above = reverse(cumsum(reverse(pdf.(dn, t)))) .* step(t)   # chance tmin >= t
        sum(pdf.(dx, t) .* above) * step(t)
    end

    # Months of all years where tmin may reach tmax on 1% of days or more
    # (year, month, expected days), most days first; and the expected days
    # over all years
    function t_crossings(t_all)
        rows, total = NamedTuple[], 0.0
        for (x, n) ∈ zip(t_all["tmax"], t_all["tmin"]), m ∈ 1:12
            p = p_cross((x.μ[m], x.σ[m], x.γ[m]), (n.μ[m], n.σ[m], n.γ[m]))
            total += p * x.mdays[m]
            p >= 0.01 && push!(rows, (year=x.year, month=m, days=p * x.mdays[m]))
        end
        sort!(rows; by=r -> -r.days), total
    end

    # Caution about days where tmin may reach tmax (nothing if none likely)
    function t_cross_note(t_all)
        rows, total = t_crossings(t_all)
        isempty(rows) && return nothing
        worst = join(["$(MONTHS[r.month]) $(r.year) ($(round(r.days; digits=1)) days)"
                      for r ∈ first(rows, 3)], ", ")
        @htl("""<p style="color: #b52f2f;"><b>Caution:</b> tmin may reach tmax on about
             $(round(total; digits=1)) day(s) over all years, in $(length(rows)) month(s) where
             it is likely on 1% of days or more; most in $(worst). MsiaGen swaps tmin and tmax
             on such days. To avoid it, widen the gap between the tmax and tmin means, or
             lower their sd, in those months.</p>""")
    end

    # Table of the resulting monthly values of year s; `monthly`: the month
    # sliders; `t_all`: every year of both variables, for the tmin/tmax check
    function t_table(s, v, monthly, t_all)
        cell(value, key, m) = "$(value) ($(fmt_change(monthly[Symbol("$(v)_$(key)$(m)")])))"
        rows = map(1:12) do m
            @htl("""<tr><td><b>$(MONTHS[m])</b></td>
                    <td>$(cell(round(s.μ[m]; digits=2), "dm", m))</td>
                    <td>$(cell(s.σ[m], "sd", m))</td>
                    <td>$(cell(s.ρ[m], "rl", m))</td>
                    <td>$(cell(s.γ[m], "sk", m))</td></tr>""")
        end
        @htl("""
        $(CENTRED)
        <p><b>Resulting monthly values of $(v), $(s.year)</b>:
           value (month slider change). Whole year: annual mean $(round(s.M; digits=2)) °C,
           amplitude $(round(s.A; digits=2)) °C, warmest day $(show_day(s.peak)).</p>
        <table class="centred">
        <tr><th>Month</th><th>mean (°C)</th><th>sd (°C)</th><th>rlag</th><th>skew</th></tr>
        $(rows)
        </table>
        $(t_cross_note(t_all))
        """)
    end
end;

# ╔═╡ feb14a42-bf45-47f8-a172-bfad503f1db9
# End-year targets of the whole-year values, and the path to them. Targets
# start at the start-year values (section 1).
@bind ttarget t_target_widget(twhole)

# ╔═╡ 8fc48e5a-3ca5-4279-bf20-b06a237784fa
# temperature statistics of every year, for tmax and tmin
t_all = all_t_stats(years, twhole, ttarget, tmonthly);

# ╔═╡ edfd7811-cdc2-4a9d-bdc3-99fad649447c
# the start year of the chosen variable, shown in the chart and the table
begin
    t_shown = first(t_all[tvar])
    t_points = t_sample_year(t_shown, t_draws)
end;

# ╔═╡ a88f601f-0b2a-4cad-8818-2054f575d03c
t_chart(t_shown, t_points, tvar,
        "$(tvar) for $(t_shown.year): illustrative days around your monthly statistics";
        other=first(t_all[other_tvar(tvar)])) |> WideCell

# ╔═╡ d591c92f-0560-410d-b8d7-f431b6e45396
t_table(t_shown, tvar, tmonthly, t_all)

# ╔═╡ 86c421d6-d3a0-421f-a6cb-8a856411a21a
# Preview of the end year, with the start year's monthly means for reference
let e = last(t_all[tvar])
    t_chart(e, t_sample_year(e, t_draws), tvar,
            "Preview: $(tvar) for $(e.year), the end year"; ref=t_shown.μ,
            other=last(t_all[other_tvar(tvar)])) |> WideCell
end

# ╔═╡ 7b58bb8a-4f4d-49e1-90dc-04b31f3529da
md"""
# 2. Wind speed

MsiaGen needs three statistics per month for wind: **mean**, **sd** and **rlag**. There is no `skew`: MsiaGen draws the daily departures from a Weibull distribution whose shape follows from `sd` relative to `mean` (a larger `sd`/`mean` gives more gusty spikes). MsiaGen never lets wind fall below **0.1 m/s**; with a low mean and a high `sd`, many days sit on that floor and the real average comes out higher than set.

Start values when no weather file is loaded (typical of Malaysian sites; the windiest time is the northeast monsoon):

| annual mean (m/s) | amplitude (m/s) | windiest day | sd (m/s) | rlag |
|---|---|---|---|---|
| 1.6 | 0.2 | 45 (mid-February) | 0.3 | 0.2 |
"""

# ╔═╡ d6ece287-169f-4eb2-acd8-c60521e0339b
begin
    # Start values of the whole year, used when no weather file is loaded:
    # medians of the 23 sites' stats files
    const W_START = (M=1.6, A=0.2, peak=45, sd=0.3, rlag=0.2)
    # valid limits of each statistic (MsiaGen's wind floor is 0.1 m/s)
    const W_LIMITS = (mean=(0.1, 10.0), sd=(0.05, 2.0), rlag=(0.0, 0.95))
    const W_FLOOR = 0.1
    # slider values: whole year, and month changes (d...)
    const W_RANGE = (M=steps(W_LIMITS.mean..., 0.05), A=steps(0.0, 3.0, 0.05), peak=1:365,
                     sd=steps(W_LIMITS.sd..., 0.01), rlag=steps(W_LIMITS.rlag..., 0.01),
                     dm=steps(-2.0, 2.0, 0.05), dsd=steps(-1.0, 1.0, 0.01),
                     drl=steps(-1, 1, 0.1))
end;

# ╔═╡ 9a83612e-c6eb-4b3b-a57f-a5c2e0252907
# Start positions of the wind sliders: from the loaded file's average year,
# or the start values with no month changes (see T_INIT)
W_INIT = let b = baseline isa NamedTuple ? baseline.wind : nothing
    if isnothing(b)
        (whole=W_START, dm=zeros(12), dsd=zeros(12), drl=zeros(12))
    else
        r = W_RANGE
        f = fit_curve(b.μ)
        M, A = snap(r.M, f.M), snap(r.A, f.A)
        sd, rl = snap(r.sd, mean(b.σ)), snap(r.rlag, mean(b.ρ))
        (whole=(M=M, A=A, peak=f.peak, sd=sd, rlag=rl),
         dm=snap.(Ref(r.dm), b.μ .- curve_months(M, A, f.peak)),
         dsd=snap.(Ref(r.dsd), b.σ .- sd), drl=snap.(Ref(r.drl), b.ρ .- rl))
    end
end;

# ╔═╡ 64e9110e-61f0-4ae0-9998-4af969a6b8f2
Markdown.parse("## a. Set start year ($(start_year))")

# ╔═╡ 5ec098d8-89b9-447b-9d26-1ccecbfd555b
@bind wwhole PlutoUI.combine() do Child
    d, r = W_INIT.whole, W_RANGE
    @htl("""
    <table>
    <tr><td>Annual mean (m/s)</td>
        <td>$(Child(:M, PlutoUI.Slider(r.M; default=d.M, show_value=true)))</td></tr>
    <tr><td>Amplitude (m/s)</td>
        <td>$(Child(:A, PlutoUI.Slider(r.A; default=d.A, show_value=true)))</td></tr>
    <tr><td>Windiest day (day of year)</td>
        <td>$(Child(:peak, PlutoUI.Slider(r.peak; default=d.peak, show_value=true)))</td></tr>
    <tr><td>sd (m/s)</td>
        <td>$(Child(:sd, PlutoUI.Slider(r.sd; default=d.sd, show_value=true)))</td></tr>
    <tr><td>rlag</td>
        <td>$(Child(:rlag, PlutoUI.Slider(r.rlag; default=d.rlag, show_value=true)))</td></tr>
    </table>
    """)
end

# ╔═╡ 488868ca-6759-430f-a07a-3283e29a696c
@bind wnew CounterButton("New sample of points")

# ╔═╡ 0b572148-a8dd-4899-9b9c-fe5734956bd8
md"""
## b. Each month
"""

# ╔═╡ 0574ff20-27cc-4ed7-947b-7f085b898452
# Each slider is a change from the whole-year value of the start year.
# Results beyond W_LIMITS are capped (and reported in the output).
@bind wmonthly PlutoUI.combine() do Child
    d, r = W_INIT, W_RANGE
    rows = map(1:12) do m
        @htl("""
        <tr>
          <td><b>$(MONTHS[m])</b></td>
          <td>$(Child(Symbol("dm$m"), PlutoUI.Slider(r.dm; default=d.dm[m], show_value=true)))</td>
          <td>$(Child(Symbol("sd$m"), PlutoUI.Slider(r.dsd; default=d.dsd[m], show_value=true)))</td>
          <td>$(Child(Symbol("rl$m"), PlutoUI.Slider(r.drl; default=d.drl[m], show_value=true)))</td>
        </tr>
        """)
    end
    @htl("""
    <table>
    <tr><th>Month</th><th>mean change (m/s)</th><th>sd change (m/s)</th><th>rlag change</th></tr>
    $(rows)
    </table>
    """)
end

# ╔═╡ 94154937-7998-4082-863e-b754759cb156
let e = max(start_year, end_year)
    Markdown.parse("""
    ## c. Set last year ($(e))

    Set the targets for $(e), and the path the values take from $(start_year) to get there; the years between follow that path, in all months. Targets start at the $(start_year) values, so changing section 1 resets them. The preview below shows $(e).
    """)
end

# ╔═╡ 561fa940-4e23-4362-9600-ff88dada2f33
begin
    # Monthly wind statistics of year y, built like temperature: whole-year
    # values moved towards the end-year targets along the chosen path, a
    # seasonal curve for the monthly means, then the month sliders. Values
    # beyond W_LIMITS are capped, and the year is flagged.
    # `w`, `t` and `monthly` are the slider values of sections 1, 3 and 2.
    function w_year_stats(y, w, t, monthly)
        mc(k) = [monthly[Symbol("$(k)$(m)")] for m ∈ 1:12]
        g = path_g(t.path, frac(y))
        M = toward(w.M, t.M, g)
        A = max(0.0, toward(w.A, t.A, g))
        peak = toward_day(w.peak, t.peak, g)
        cal = calendar(y)
        curve = M .+ A .* cos.(2π .* ((1:cal.ndays) .- peak) ./ cal.ndays)
        raw = (mean=[mean(curve[cal.month_of .== m]) for m ∈ 1:12] .+ mc("dm"),
               sd=toward(w.sd, t.sd, g) .+ mc("sd"),
               rlag=toward(w.rlag, t.rlag, g) .+ mc("rl"))
        capped(k) = round.(clamp.(raw[k], W_LIMITS[k]...); digits=3)
        (year=y, M=M, A=A, peak=peak,
         μ=capped(:mean), σ=capped(:sd), ρ=capped(:rlag),
         mdays=cal.mdays, ndays=cal.ndays, month_of=cal.month_of,
         clamped=any(any(raw[k] .!= clamp.(raw[k], W_LIMITS[k]...)) for k ∈ keys(raw)))
    end

    # One illustrative year of wind: MsiaGen's wind model, month by month
    # (genwind.jl: wind_dist!), with fixed random numbers. Departures come from
    # a Weibull distribution with the month's mean and the SD left after the
    # lag-1 term; wind never falls below 0.1 m/s.
    function w_sample_year(s, u)
        x = zeros(s.ndays)
        t0 = 1
        for m ∈ 1:12
            idx = t0:t0+s.mdays[m]-1
            sde = s.σ[m] * sqrt(1 - s.ρ[m]^2)
            shape = (sde / s.μ[m])^-1.086
            scale = s.μ[m] / gamma(1 + 1 / shape)
            e = quantile.(Weibull(shape, scale), u[idx]) .- s.μ[m]
            c = s.μ[m] * (1 - s.ρ[m])
            for (k, t) ∈ enumerate(idx)
                prev = (t == 1) ? s.μ[m] : x[t-1]
                x[t] = max(W_FLOOR, c + s.ρ[m] * prev + e[k])
            end
            t0 += s.mdays[m]
        end
        x
    end

    # Chart of one year's illustrative wind around its monthly statistics;
    # `ref`: monthly means of another year, drawn as a dashed reference
    w_chart(s, pts, title; ref=nothing) =
        day_chart(s, pts, title, "wind speed (m/s)", "m/s"; ref=ref, floor_at=W_FLOOR,
                  floor_label="0.1 m/s floor")

    all_w_stats(years, whole, target, monthly) = [w_year_stats(y, whole, target, monthly) for y ∈ years]

    # Sliders of the end-year targets of the whole-year values and the path
    # to them (section 3), starting at the start-year values `whole`
    function w_target_widget(whole)
            PlutoUI.combine() do Child
            r = W_RANGE
            now(k) = whole[k]
            sl(k) = Child(k, PlutoUI.Slider(r[k]; default=now(k), show_value=true))
            row(label, k) = @htl("""<tr><td>$(label)</td><td>$(sl(k))</td><td>$(now(k))</td></tr>""")
            @htl("""
            <p><b>Path:</b> $(Child(:path, Select(PATHS)))</p>
            <table>
            <tr><th></th><th>Target in end year</th><th>Start year</th></tr>
            $(row("Annual mean (m/s)", :M))
            $(row("Amplitude (m/s)", :A))
            $(row("Windiest day (day of year)", :peak))
            </table>
            <details><summary>More settings: sd, rlag</summary>
            <table>
            <tr><th></th><th>Target in end year</th><th>Start year</th></tr>
            $(row("sd (m/s)", :sd))
            $(row("rlag", :rlag))
            </table>
            </details>
            """)
        end
    end

    # Fixed random numbers for the floor check, so it does not flicker
    const W_CHECK_U = rand(MersenneTwister(99), 3000)

    # Share of days on the 0.1 m/s floor in a month with wind mean μ, sd σ
    # and rlag ρ: MsiaGen's recursion run over 3000 days. MsiaGen's search
    # for the best month adds some more (up to about 30% more, in tests).
    function w_floor_share(μ, σ, ρ)
        sde = σ * sqrt(1 - ρ^2)
        shape = (sde / μ)^-1.086
        scale = μ / gamma(1 + 1 / shape)
        c, x, n = μ * (1 - ρ), μ, 0
        for u ∈ W_CHECK_U
            x = max(W_FLOOR, c + ρ * x + scale * (-log1p(-u))^(1 / shape) - μ)
            n += x <= W_FLOOR
        end
        n / length(W_CHECK_U)
    end

    # Months of all years (`w_all`) where wind may sit on the floor on 1% of
    # days or more: year, month and share of days, most first
    w_floored(w_all) = sort!([(year=s.year, month=m, share=f) for s ∈ w_all for m ∈ 1:12
                              for f ∈ (w_floor_share(s.μ[m], s.σ[m], s.ρ[m]),) if f >= 0.01];
                             by=r -> -r.share)

    # Caution about wind on the floor (nothing if not likely)
    function w_floor_note(w_all)
        rows = w_floored(w_all)
        isempty(rows) && return nothing
        worst = join(["$(MONTHS[r.month]) $(r.year) ($(round(Int, 100 * r.share))%)"
                      for r ∈ first(rows, 3)], ", ")
        @htl("""<p style="color: #b52f2f;"><b>Caution:</b> in $(length(rows)) month(s), wind
             may fall to the 0.1 m/s floor on at least 1% of days; most in $(worst) (at least
             that share of days). There, the generated wind sits on the floor more often than
             real wind does, and its mean comes out higher than set. To avoid it, lower the sd
             or raise the mean in those months.</p>""")
    end

    # Table of the resulting monthly values of year s; `monthly`: the month
    # sliders; `w_all`: every year, for the floor check
    function w_table(s, monthly, w_all)
        cell(value, key, m) = "$(value) ($(fmt_change(monthly[Symbol("$(key)$(m)")])))"
        rows = map(1:12) do m
            @htl("""<tr><td><b>$(MONTHS[m])</b></td>
                    <td>$(cell(round(s.μ[m]; digits=2), "dm", m))</td>
                    <td>$(cell(s.σ[m], "sd", m))</td>
                    <td>$(cell(s.ρ[m], "rl", m))</td></tr>""")
        end
        note = w_floor_note(w_all)
        @htl("""
        $(CENTRED)
        <p><b>Resulting monthly values of wind speed, $(s.year)</b>:
           value (month slider change). Whole year: annual mean $(round(s.M; digits=2)) m/s,
           amplitude $(round(s.A; digits=2)) m/s, windiest day $(show_day(s.peak)).</p>
        <table class="centred">
        <tr><th>Month</th><th>mean (m/s)</th><th>sd (m/s)</th><th>rlag</th></tr>
        $(rows)
        </table>
        $(note)
        """)
    end
end;

# ╔═╡ 1f8a6cf2-8e88-4ec3-89e0-357e8e76112f
# End-year targets of the whole-year values, and the path to them. Targets
# start at the start-year values (section 1).
@bind wtarget w_target_widget(wwhole)

# ╔═╡ 1ed9263f-1d99-4105-beb2-5a5a6eaa477f
w_all = all_w_stats(years, wwhole, wtarget, wmonthly);

# ╔═╡ cb235b2e-649e-4679-a53f-dce0ed1a4b0a
# fixed random numbers, so points move smoothly when a slider changes
w_draws = let rng = MersenneTwister(1234 + wnew)
    rand(rng, 366)
end;

# ╔═╡ 0ac15b87-47c7-4b4d-8d79-54d9a9361779
# the start year, shown in the chart and the table
begin
    w_shown = first(w_all)
    w_points = w_sample_year(w_shown, w_draws)
end;

# ╔═╡ 75d39779-4ffa-4128-8b11-7ec6ee171f1f
w_chart(w_shown, w_points,
        "Wind speed for $(w_shown.year): illustrative days around your monthly statistics") |> WideCell

# ╔═╡ 2b321dcb-93b2-421d-a75c-747999faa12b
w_table(w_shown, wmonthly, w_all)

# ╔═╡ 4f563fdd-9a55-4f0c-a86e-90135154a104
# Preview of the end year, with the start year's monthly means for reference
let e = last(w_all)
    w_chart(e, w_sample_year(e, w_draws), "Preview: wind speed for $(e.year), the end year";
            ref=w_shown.μ) |> WideCell
end

# ╔═╡ 6bbc7b6c-e62d-4244-82e1-422cb6c002c3
md"""
# 3. Rainfall

MsiaGen needs three statistics per month for rain: the **total rainfall**, **pww** and **pwd**:

- Wet days: Higher `pww` = longer wet spells.
- Dry days: Lower `pwd` = longer dry spells.

Together they set how many days are wet, and how wet and dry days bunch together. The designer shows their effect as **expected rain days** and the **average wet and dry spell** lengths.

Start values when no weather file is loaded:

| annual rainfall (mm) | seasonality (%) | driest day | pww | pwd | rain days in a 30-day month |
|---|---|---|---|---|---|
| 2460 | 30 | 183 | 0.58 | 0.42 | 15 |

**Seasonality** is how much wetter the wettest time of year is than average, and drier the driest: at 30%, the driest months get about 0.7× and the wettest about 1.3× the average monthly rainfall. **Driest day** is the middle of the dry season. Mid-year suits the west coast; in the south and on the east coast the driest time is usually February–March (day 45–75).

Expected rain days follow MsiaGen exactly: days in the month × `pwd` / (1 − `pww` + `pwd`), rounded to the nearest day (at least 1 in a month with rain). The calendar strip is **illustrative**: wet and dry days drawn with MsiaGen's own rain model so you can see what the numbers mean. It is not the generated weather.
"""

# ╔═╡ c322bc84-d17d-460d-8694-ff964e96f44b
begin
    # Start values of the whole year, used when no weather file is loaded.
    # Annual rainfall, pww and pwd are medians of the 23 sites' stats files;
    # seasonality (%) and the driest day (day of year) set the seasonal curve.
    const R_START = (total=2460.0, seas=30.0, dry=183, pww=0.58, pwd=0.42)
    # valid limits: rainfall 0 or more; seasonality 0-90 % (so no month gets
    # negative rain); pww and pwd strictly between 0 and 1
    const R_LIMITS = (total=(0.0, Inf), seas=(0.0, 90.0), pww=(0.01, 0.99), pwd=(0.01, 0.99))
    # slider values: whole year, and month changes (d...)
    const R_RANGE = (total=steps(100.0, 6000.0, 10.0), seas=steps(R_LIMITS.seas..., 5.0),
                     dry=1:365, pww=steps(R_LIMITS.pww..., 0.01), pwd=steps(R_LIMITS.pwd..., 0.01),
                     dr=steps(-100.0, 200.0, 1.0), dww=steps(-0.5, 0.5, 0.01),
                     dwd=steps(-0.5, 0.5, 0.01))
    # daily rain classes (mm) for the calendar strip
    const RAIN_BOUNDS = [0.0, 5.0, 10.0, 20.0, 50.0]
    const RAIN_COLORS = ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
end;

# ╔═╡ 5b47c4fa-56b2-47f3-ad93-c63da0ff47b2
# Start positions of the rain sliders: from the loaded file's average year
# (annual total; seasonal curve fitted to its monthly totals; each month's
# rainfall change in % holds the remainder), or the start values with no
# month changes
R_INIT = let b = baseline isa NamedTuple ? baseline.rain : nothing
    if isnothing(b)
        (whole=R_START, dr=zeros(12), dww=zeros(12), dwd=zeros(12))
    else
        r = R_RANGE
        f = fit_rain(b.total)
        total, seas = snap(r.total, sum(b.total)), snap(r.seas, f.seas)
        pww, pwd = snap(r.pww, mean(b.pww)), snap(r.pwd, mean(b.pwd))
        expected = total .* rain_shares(seas, f.dry)
        (whole=(total=total, seas=seas, dry=f.dry, pww=pww, pwd=pwd),
         dr=snap.(Ref(r.dr), 100 .* (b.total ./ expected .- 1)),
         dww=snap.(Ref(r.dww), b.pww .- pww), dwd=snap.(Ref(r.dwd), b.pwd .- pwd))
    end
end;

# ╔═╡ 43e129a8-1208-4b4a-8669-ace699188310
Markdown.parse("## a. Set start year ($(start_year))")

# ╔═╡ c64cc811-9f2e-4529-ac78-be3fbfcdbad2
@bind rwhole PlutoUI.combine() do Child
    d, r = R_INIT.whole, R_RANGE
    @htl("""
    <table>
    <tr><td>Annual rainfall (mm)</td>
        <td>$(Child(:total, PlutoUI.Slider(r.total; default=d.total, show_value=true)))</td></tr>
    <tr><td>Seasonality (%)</td>
        <td>$(Child(:seas, PlutoUI.Slider(r.seas; default=d.seas, show_value=true)))</td></tr>
    <tr><td>Driest day (day of year)</td>
        <td>$(Child(:dry, PlutoUI.Slider(r.dry; default=d.dry, show_value=true)))</td></tr>
    <tr><td>pww</td>
        <td>$(Child(:pww, PlutoUI.Slider(r.pww; default=d.pww, show_value=true)))</td></tr>
    <tr><td>pwd</td>
        <td>$(Child(:pwd, PlutoUI.Slider(r.pwd; default=d.pwd, show_value=true)))</td></tr>
    </table>
    """)
end

# ╔═╡ 14157358-7ff1-4a10-ac6f-71a3258ea780
@bind rnew CounterButton("Re-sample rain days")

# ╔═╡ d4682d65-c35a-44e1-b6a3-496e3b3cbb2f
md"""
## b. Each month
"""

# ╔═╡ 96299f6a-86d7-4c94-8e6c-6b7ebe001953
# Each slider is a change from the whole-year value of the start year:
# rainfall in % of the month's share, pww and pwd as additions.
# Results beyond R_LIMITS are capped (and reported in the output).
@bind rmonthly PlutoUI.combine() do Child
    d, r = R_INIT, R_RANGE
    rows = map(1:12) do m
        @htl("""
        <tr>
          <td><b>$(MONTHS[m])</b></td>
          <td>$(Child(Symbol("dr$m"), PlutoUI.Slider(r.dr; default=d.dr[m], show_value=true)))</td>
          <td>$(Child(Symbol("ww$m"), PlutoUI.Slider(r.dww; default=d.dww[m], show_value=true)))</td>
          <td>$(Child(Symbol("wd$m"), PlutoUI.Slider(r.dwd; default=d.dwd[m], show_value=true)))</td>
        </tr>
        """)
    end
    @htl("""
    <table>
    <tr><th>Month</th><th>rainfall change (%)</th><th>pww change</th><th>pwd change</th></tr>
    $(rows)
    </table>
    """)
end

# ╔═╡ ce6ecf90-981c-4037-91fd-66bbcdec3178
let e = max(start_year, end_year)
    Markdown.parse("""
    ## c. Set last year ($(e))

    Set the targets for $(e), and the path the values take from $(start_year) to get there; the years between follow that path, in all months. Targets start at the $(start_year) values, so changing section 1 resets them. The preview below shows $(e).
    """)
end

# ╔═╡ c538e51c-8365-4070-b471-b1d4c49a3511
begin
    # Expected rain days in a month, as MsiaGen sets them (genrain.jl:
    # gen_rain_month): days × pwd / (1 - pww + pwd), rounded, at least 1 if
    # the month has rain and at most the days in the month. Average spell
    # lengths of the wet/dry chain.
    begin
        wet_fraction(pww, pwd) = (d = 1 - pww + pwd; d ≈ 0 ? 1.0 : pwd / d)
        rain_days(days, pww, pwd, total) =
            clamp(round(Int, days * wet_fraction(pww, pwd)), total > 0 ? 1 : 0, days)
        wet_spell(pww) = 1 / (1 - pww)
        dry_spell(pwd) = 1 / pwd
    end

    # Monthly rain statistics of year y. The whole-year values move from the
    # start year (section 1) towards the end-year targets (section 3) along
    # the chosen path. The annual rainfall is spread over the months by the
    # year's seasonal curve, lowest on the driest day; the month sliders are
    # then applied. Values beyond R_LIMITS are capped, and the year is flagged.
    # `w`, `t` and `monthly` are the slider values of sections 1, 3 and 2.
    function r_year_stats(y, w, t, monthly)
        mc(k) = [monthly[Symbol("$(k)$(m)")] for m ∈ 1:12]
        g = path_g(t.path, frac(y))
        annual = toward(w.total, t.total, g)
        seas = clamp(toward(w.seas, t.seas, g), R_LIMITS.seas...)
        dry = toward_day(w.dry, t.dry, g)
        cal = calendar(y)
        curve = 1 .- (seas / 100) .* cos.(2π .* ((1:cal.ndays) .- dry) ./ cal.ndays)
        share = [sum(curve[cal.month_of .== m]) for m ∈ 1:12] ./ sum(curve)
        raw = (total=annual .* share .* (1 .+ mc("dr") ./ 100),
               pww=toward(w.pww, t.pww, g) .+ mc("ww"),
               pwd=toward(w.pwd, t.pwd, g) .+ mc("wd"))
        capped(k) = round.(clamp.(raw[k], R_LIMITS[k]...); digits=3)
        total, pww, pwd = capped(:total), capped(:pww), capped(:pwd)
        (year=y, seas=seas, dry=dry, total=total, pww=pww, pwd=pwd, mdays=cal.mdays,
         raindays=rain_days.(cal.mdays, pww, pwd, total),
         wetspell=wet_spell.(pww), dryspell=dry_spell.(pwd),
         clamped=any(any(raw[k] .!= clamp.(raw[k], R_LIMITS[k]...)) for k ∈ keys(raw)))
    end    

    # Illustrative daily rain of each month of year s, following MsiaGen
    # (genrain.jl): n = expected rain days; n amounts from a Gamma
    # distribution whose shape comes from MsiaGen's GEV; wet days placed by
    # the pww/pwd chain (extra chain-wet days stay dry, missing ones are taken
    # from dry days). Amounts are scaled so each month's total is exactly the
    # set total.
    function make_strip(s, draws)
        gev = GeneralizedExtremeValue(0.50, 0.17, 0.14)
        wet_prev = false
        map(1:12) do m
            days = s.mdays[m]
            n = s.raindays[m]
            rain = zeros(days)
            (n == 0 || s.total[m] <= 0) && return rain
            wet = falses(days)
            prev = m == 1 ? draws.chain[m, 1] <= wet_fraction(s.pww[m], s.pwd[m]) : wet_prev
            for d ∈ 1:days
                wet[d] = draws.chain[m, d] <= (prev ? s.pww[m] : s.pwd[m])
                prev = wet[d]
            end
            pos = findall(wet)
            if length(pos) >= n
                pos = pos[1:n]
            else
                dry = findall(.!wet)
                dry = dry[sortperm(draws.extra[m, dry])]
                pos = sort([pos; dry[1:n-length(pos)]])
            end
            k = max(0.05, quantile(gev, clamp(draws.shape[m], 0.05, 0.95)))
            amt = quantile.(Gamma(k, 1.0), clamp.(draws.amount[m, 1:n], 1e-6, 1 - 1e-6))
            rain[pos] = amt .* (s.total[m] / sum(amt))
            wet_prev = rain[end] > 0
            rain
        end
    end

    # Illustrative daily rain (left) and monthly rainfall (right) of year s,
    # one row per month; `ref`: monthly totals of another year, marked on the
    # bars as a reference. Hover a day to see its rain.
    function r_chart(s, strip, title; ref=nothing)
        W, H = 1100, 500
        T, B = 44, 88                       # top and bottom margins
        gx, gw = 50, 640                    # day grid: left edge, width
        bx, bw = 760, 230                   # bars: left edge, width
        rowh, cw = (H - T - B) / 12, gw / 31
        top = maximum(s.total)
        isnothing(ref) || (top = max(top, maximum(ref)))
        xmax = 1.6 * max(1.0, top)
        BX(v) = bx + v / xmax * bw
        io = IOBuffer()
        print(io, """<text x="$gx" y="26" font-size="17" font-weight="600" fill="$INK">$(svg_text(title))</text>""",
                  """<text x="$bx" y="26" font-size="15" font-weight="600" fill="$INK">Monthly rainfall</text>""")
        for v ∈ nice_ticks(0, xmax; n=3)
            x = r1(BX(v))
            print(io, """<line x1="$x" x2="$x" y1="$T" y2="$(H - B)" stroke="$GRID"/>""",
                      """<text x="$x" y="$(H - B + 22)" font-size="14" fill="$INK2" text-anchor="middle">$(tick_label(v))</text>""")
        end
        print(io, """<line x1="$bx" x2="$bx" y1="$T" y2="$(H - B)" stroke="$GRID"/>""",
                  """<text x="$(bx + bw / 2)" y="$(H - B + 44)" font-size="14" fill="$INK2" text-anchor="middle">rainfall (mm)</text>""")
        for m ∈ 1:12
            y = T + (m - 1) * rowh
            print(io, """<text x="$(gx - 8)" y="$(r1(y + rowh / 2 + 5))" font-size="14" fill="$INK2" text-anchor="end">$(MONTHS[m])</text>""")
            for d ∈ 1:s.mdays[m]
                r = strip[m][d]
                k = r > 0 ? searchsortedlast(RAIN_BOUNDS, r) : 0
                tip = r > 0 ? "$(round(r; digits=1)) mm" : "dry"
                print(io, """<rect x="$(r1(gx + (d - 1) * cw + 1))" y="$(r1(y + 1))" width="$(r1(cw - 2))" height="$(r1(rowh - 2))" fill="$(RAIN_COLORS[k + 1])"><title>$(MONTHS[m]) $d: $tip</title></rect>""")
            end
            print(io, """<rect x="$bx" y="$(r1(y + 0.15rowh))" width="$(r1(BX(s.total[m]) - bx))" height="$(r1(0.7rowh))" fill="$BLUE"/>""",
                      """<text x="$(r1(BX(s.total[m]) + 6))" y="$(r1(y + rowh / 2 + 5))" font-size="14" fill="$INK2">$(round(Int, s.total[m])) mm, $(s.raindays[m]) d</text>""")
            if !isnothing(ref)
                x = r1(BX(ref[m]))
                print(io, """<line x1="$x" x2="$x" y1="$(r1(y + 0.1rowh))" y2="$(r1(y + 0.9rowh))" stroke="$INK2" stroke-width="2.5"/>""")
            end
        end
        for d ∈ [1; 5:5:30]
            print(io, """<text x="$(r1(gx + (d - 0.5) * cw))" y="$(H - B + 22)" font-size="14" fill="$INK2" text-anchor="middle">$d</text>""")
        end
        print(io, """<text x="$(gx + gw / 2)" y="$(H - B + 44)" font-size="14" fill="$INK2" text-anchor="middle">day of month</text>""")
        svg_legend(io, gx, H - 18, [(:box, c, l) for (c, l) ∈
                   zip(RAIN_COLORS, ["dry", "< 5", "5–10", "10–20", "20–50", "≥ 50"])]; title="mm/day")
        isnothing(ref) || svg_legend(io, bx, H - 18, [(:line, INK2, "start-year total")])
        svg(String(take!(io)), W, H)
    end

    all_r_stats(years, whole, target, monthly) = [r_year_stats(y, whole, target, monthly) for y ∈ years]

    # Sliders of the end-year targets of the whole-year values and the path
    # to them (section 3), starting at the start-year values `whole`
    function r_target_widget(whole)
            PlutoUI.combine() do Child
            r = R_RANGE
            now(k) = whole[k]
            sl(k) = Child(k, PlutoUI.Slider(r[k]; default=now(k), show_value=true))
            row(label, k) = @htl("""<tr><td>$(label)</td><td>$(sl(k))</td><td>$(now(k))</td></tr>""")
            @htl("""
            <p><b>Path:</b> $(Child(:path, Select(PATHS)))</p>
            <table>
            <tr><th></th><th>Target in end year</th><th>Start year</th></tr>
            $(row("Annual rainfall (mm)", :total))
            $(row("Seasonality (%)", :seas))
            $(row("Driest day (day of year)", :dry))
            </table>
            <details><summary>More settings: pww, pwd</summary>
            <table>
            <tr><th></th><th>Target in end year</th><th>Start year</th></tr>
            $(row("pww", :pww))
            $(row("pwd", :pwd))
            </table>
            </details>
            """)
        end
    end

    # Table of the resulting monthly values of year s; `monthly`: the month
    # sliders
    function r_table(s, monthly)
        ch(key, m) = monthly[Symbol("$(key)$(m)")]
        rows = map(1:12) do m
            @htl("""<tr><td><b>$(MONTHS[m])</b></td>
                    <td>$(round(Int, s.total[m])) ($(fmt_change(ch("dr", m)))%)</td>
                    <td>$(s.pww[m]) ($(fmt_change(ch("ww", m))))</td>
                    <td>$(s.pwd[m]) ($(fmt_change(ch("wd", m))))</td>
                    <td>$(s.raindays[m])</td>
                    <td>$(round(s.wetspell[m]; digits=1))</td>
                    <td>$(round(s.dryspell[m]; digits=1))</td></tr>""")
        end
        @htl("""
        $(CENTRED)
        <p><b>Resulting monthly values of rainfall, $(s.year)</b>: value (month slider change).
           Whole year: annual rainfall $(round(Int, sum(s.total))) mm,
           seasonality $(round(Int, s.seas))%, driest day $(show_day(s.dry)).</p>
        <table class="centred">
        <tr><th>Month</th><th>rainfall (mm)</th><th>pww</th><th>pwd</th>
            <th>rain days</th><th>avg wet spell (days)</th><th>avg dry spell (days)</th></tr>
        $(rows)
        </table>
        """)
    end
end;

# ╔═╡ ad4f0f49-552f-4e2a-9ccc-2ae81c2e6584
# End-year targets of the whole-year values, and the path to them. Targets
# start at the start-year values (section 1).
@bind rtarget r_target_widget(rwhole)

# ╔═╡ 2146995a-4493-4b51-8458-a266995a566e
begin
    r_all = all_r_stats(years, rwhole, rtarget, rmonthly)
    r_shown = first(r_all)     # the start year, shown in the chart and table
end;

# ╔═╡ a4d95768-97e4-4f28-9244-26918a4ba23f
r_table(r_shown, rmonthly)

# ╔═╡ cceb8405-6585-4798-b665-a90fb032b1af
# fixed random numbers, so the strip changes smoothly when a slider moves
r_draws = let rng = MersenneTwister(1234 + rnew)
    (chain=rand(rng, 12, 31),     # wet/dry chain
     extra=rand(rng, 12, 31),     # order of dry days turned wet, if needed
     amount=rand(rng, 12, 31),    # rain amounts
     shape=rand(rng, 12))         # Gamma shape of each month
end;

# ╔═╡ 73675405-fcb9-46ae-b4a3-a85a0a5da77d
# Preview of the end year, with the start year's monthly totals for reference
let e = last(r_all)
    r_chart(e, make_strip(e, r_draws), "Preview: rain days for $(e.year), the end year";
            ref=r_shown.total) |> WideCell
end

# ╔═╡ 3d93cc35-b1eb-4efd-b7ba-4279c171692f
# illustrative daily rain of the start year
r_strip = make_strip(r_shown, r_draws);

# ╔═╡ 76dc5247-761a-4391-9b85-e49d6877f40c
r_chart(r_shown, r_strip, "Illustrative rain days for $(r_shown.year)") |> WideCell

# ╔═╡ b7ab44db-4cdb-4b75-a322-c8e997dcbc65
md"""
# 4. Output for MsiaGen

$(@bind ok CounterButton("Generate"))
"""

# ╔═╡ c071bc59-7701-4533-9e3f-76f451b2f140
# One MsiaGen stats row per year: tmin, tmax, wind, then rain (the order
# MsiaGen itself writes). Annual (0) columns: temperature and wind mean are
# day-weighted means of the months, rain total is the sum; the other annual
# values are averages of the 12 months (MsiaGen uses them only for its
# progress messages and the first day of the year).
begin
    t_cols(s) = [sum(s.μ .* s.mdays) / s.ndays; s.μ;
                 mean(s.σ); s.σ; mean(s.ρ); s.ρ; mean(s.γ); s.γ]
    w_cols(s) = [sum(s.μ .* s.mdays) / s.ndays; s.μ;
                 mean(s.σ); s.σ; mean(s.ρ); s.ρ]
    r_cols(s) = [sum(s.total); s.total; mean(s.pww); s.pww; mean(s.pwd); s.pwd]
    stats_header = ["year";
        ["$(p)_$(v)$(i)" for v ∈ ("tmin", "tmax") for p ∈ ("mean", "sd", "rlag", "skew") for i ∈ 0:12];
        ["$(p)_wind$(i)" for p ∈ ("mean", "sd", "rlag") for i ∈ 0:12];
        ["$(p)$(i)" for p ∈ ("totrain", "pww", "pwd") for i ∈ 0:12]]

    # The stats file: latitude line, header, then one row per year
    function stats_csv(lat, years, t_all, w_all, r_all)
        rows = [[y; t_cols(t_all["tmin"][i]); t_cols(t_all["tmax"][i]);
                 w_cols(w_all[i]); r_cols(r_all[i])] for (i, y) ∈ enumerate(years)]
        "$(lat)\n" * join(stats_header, ",") * "\n" *
            join([join([string(Int(r[1])); string.(round.(r[2:end]; digits=4))], ",")
                  for r ∈ rows], "\n") * "\n"
    end

    # Checks: years where values were capped, and the cautions about tmin
    # reaching tmax and wind on its floor (nothing when there are none)
    function value_checks(years, t_all, w_all, r_all)
        capped_years = (temperature=sort(unique([s.year for v ∈ TVARS for s ∈ t_all[v] if s.clamped])),
                        wind=[s.year for s ∈ w_all if s.clamped],
                        rain=[s.year for s ∈ r_all if s.clamped])
        capped_years, filter(!isnothing, [t_cross_note(t_all), w_floor_note(w_all)])
    end

    # Annual summary, notes from the checks, and the stats file to download
    function results_html(site, years, t_all, w_all, r_all, capped_years, cautions,
                          csv_text, csv_name)
        notes = []
        for (part, ys) ∈ pairs(capped_years)
            isempty(ys) || push!(notes, @htl("""<p style="color: #b52f2f;"><b>Note:</b> some
                $(part) values went past their limits in $(length(ys)) year(s), $(first(ys))–$(last(ys)).
                They were capped at the limit.</p>"""))
        end
        append!(notes, cautions)
        rows = map(enumerate(years)) do (i, y)
            tx, tn, w, r = t_all["tmax"][i], t_all["tmin"][i], w_all[i], r_all[i]
            @htl("""<tr><td><b>$(y)</b></td>
                    <td>$(round(sum(tx.μ .* tx.mdays) / tx.ndays; digits=2))</td>
                    <td>$(round(sum(tn.μ .* tn.mdays) / tn.ndays; digits=2))</td>
                    <td>$(round(sum(w.μ .* w.mdays) / w.ndays; digits=2))</td>
                    <td>$(round(Int, sum(r.total)))</td>
                    <td>$(sum(r.raindays))</td></tr>""")
        end
        @htl("""
        $(CENTRED)
        <p>Annual summary of <b>$(strip(site))</b>, $(first(years))–$(last(years)) ($(length(years)) years).
           Monthly values for every year are in the CSV below.</p>
        $(notes)
        <div style="max-height: 320px; overflow-y: auto;">
        <table class="centred">
        <tr><th>Year</th><th>tmax mean (°C)</th><th>tmin mean (°C)</th><th>wind mean (m/s)</th>
            <th>rainfall (mm)</th><th>rain days</th></tr>
        $(rows)
        </table>
        </div>
        <p>In MsiaGen's stats format, to be saved as <code>$(csv_name)</code> in the site's data folder:</p>
        <pre style="white-space: pre; overflow: auto; max-height: 240px; font-size: 11px;">$(csv_text)</pre>
        $(DownloadButton(csv_text, csv_name))
        """)
    end
end;

# ╔═╡ ea084619-514a-4fa6-b7ff-ed749b27cccb
csv_text = stats_csv(lat, years, t_all, w_all, r_all);

# ╔═╡ 24156aa7-f981-412e-b05e-52e6eec5c9db
# download file name: <site>-stats.csv, the name MsiaGen reads (spaces in
# the site name become hyphens)
csv_name = "$(replace(strip(site), r"\s+" => "-"))-stats.csv";

# ╔═╡ f9ca3c6d-59a1-41cc-a955-2a8a39d4a6cc
# Checks: years where values were capped, and cautions about tmin reaching
# tmax and wind on its floor
capped_years, cautions = value_checks(years, t_all, w_all, r_all);

# ╔═╡ b7503956-1412-4e7d-855b-480d1209414f
if start_year > end_year
    md"""
    !!! warning "Start year > End year"
        Start year is after End year. Please correct the range. 
    """
elseif isempty(strip(site))
    md"""
    !!! warning "Site name not given"
        Enter the site name.
    """
elseif ok == 0
    md"_Click **Generate** to show the statistics._"
else
    results_html(site, years, t_all, w_all, r_all, capped_years, cautions, csv_text, csv_name)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Distributions = "~0.25.131"
HypertextLiteral = "~0.9.5"
PlutoUI = "~0.7.73"
SpecialFunctions = "~2.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "1f9c8ad678fde0009303c3d0c0c71980a21ab910"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "61761f58648aa7217445f24f841839b78c712232"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.3"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a958ab3a40c755563f5e1405c0846cb0446bf19d"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.131"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "086b5fbd032baf544678cc15b52b015e4c4aceb8"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.1"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "31bb6c92405c084617facc1d7ed9eb6c402d061e"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.30"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "b85e2797b2409570e84c4de46238c0ed5f6476ae"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.2"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "f9b03759e9ef463718934fbed30820b39997511e"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "ba0dc8a8a67cacac4842631f960c046e4e563675"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "5005266de4bfe50e53ff44a5cb5c540b6e47a254"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "971f04b3780c0da4edf230573fc1aadc2c6e649e"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "429071b23f4c9a13fb6582f807cc2ef454082408"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╟─695aa061-f11d-4895-a66b-093f79ebdb2e
# ╟─a2a81510-a24d-481d-a7de-a5bccbc77218
# ╟─20609a8a-99cc-4411-b1fa-8f1e2a709609
# ╟─6ad69763-051b-4d69-a99f-69c495577cab
# ╟─86792682-a96d-4ab2-a321-4199e2e5c0b1
# ╟─c6cfd927-288b-4ada-b19c-5dc6fccc2a1c
# ╟─1db20637-6907-44cc-8728-400469155d76
# ╟─1f5fb93a-a616-4cd2-8eb4-2d68e9f6a07e
# ╟─d49f9ae6-9b07-484d-a105-13b3fbac522f
# ╟─a4424332-1266-4ad7-8559-e0ba6ceb7a54
# ╟─2ac2e9cc-c555-46bc-9bf0-47bc6340192d
# ╟─cf725ce2-bebe-4b48-8405-d369e1cef355
# ╟─7a6dc0d6-5d28-4ea3-a728-3a873958d341
# ╟─a28b530c-ac28-4989-92aa-448efae7e9ad
# ╟─6d7815aa-155b-438f-999d-f1b78b356d3a
# ╟─4d2a0300-e697-4deb-8580-11e978893162
# ╟─05e5deb2-0c98-415d-a129-b8db81c24db4
# ╟─a88f601f-0b2a-4cad-8818-2054f575d03c
# ╟─192267b9-fadc-4bc5-9bb0-afa0712d1e49
# ╟─c0544ade-de13-48eb-937d-ef45ce85700c
# ╟─77485903-cd54-41ae-ad0f-460ecb1e6058
# ╟─d591c92f-0560-410d-b8d7-f431b6e45396
# ╟─3cc149d9-614e-4ad1-a2e8-fe57e395a2aa
# ╟─feb14a42-bf45-47f8-a172-bfad503f1db9
# ╟─86c421d6-d3a0-421f-a6cb-8a856411a21a
# ╟─8fc48e5a-3ca5-4279-bf20-b06a237784fa
# ╟─4b6ac062-0db7-4d6d-8649-973e913c59d8
# ╟─eea66297-900f-48a8-af1d-4881cb5fd1cd
# ╟─edfd7811-cdc2-4a9d-bdc3-99fad649447c
# ╟─7b58bb8a-4f4d-49e1-90dc-04b31f3529da
# ╟─d6ece287-169f-4eb2-acd8-c60521e0339b
# ╟─9a83612e-c6eb-4b3b-a57f-a5c2e0252907
# ╟─64e9110e-61f0-4ae0-9998-4af969a6b8f2
# ╟─5ec098d8-89b9-447b-9d26-1ccecbfd555b
# ╟─75d39779-4ffa-4128-8b11-7ec6ee171f1f
# ╟─488868ca-6759-430f-a07a-3283e29a696c
# ╟─0b572148-a8dd-4899-9b9c-fe5734956bd8
# ╟─0574ff20-27cc-4ed7-947b-7f085b898452
# ╟─2b321dcb-93b2-421d-a75c-747999faa12b
# ╟─94154937-7998-4082-863e-b754759cb156
# ╟─1f8a6cf2-8e88-4ec3-89e0-357e8e76112f
# ╟─4f563fdd-9a55-4f0c-a86e-90135154a104
# ╟─561fa940-4e23-4362-9600-ff88dada2f33
# ╟─1ed9263f-1d99-4105-beb2-5a5a6eaa477f
# ╟─cb235b2e-649e-4679-a53f-dce0ed1a4b0a
# ╟─0ac15b87-47c7-4b4d-8d79-54d9a9361779
# ╟─6bbc7b6c-e62d-4244-82e1-422cb6c002c3
# ╟─c322bc84-d17d-460d-8694-ff964e96f44b
# ╟─5b47c4fa-56b2-47f3-ad93-c63da0ff47b2
# ╟─43e129a8-1208-4b4a-8669-ace699188310
# ╟─c64cc811-9f2e-4529-ac78-be3fbfcdbad2
# ╟─76dc5247-761a-4391-9b85-e49d6877f40c
# ╟─14157358-7ff1-4a10-ac6f-71a3258ea780
# ╟─d4682d65-c35a-44e1-b6a3-496e3b3cbb2f
# ╟─96299f6a-86d7-4c94-8e6c-6b7ebe001953
# ╟─a4d95768-97e4-4f28-9244-26918a4ba23f
# ╟─ce6ecf90-981c-4037-91fd-66bbcdec3178
# ╟─ad4f0f49-552f-4e2a-9ccc-2ae81c2e6584
# ╟─73675405-fcb9-46ae-b4a3-a85a0a5da77d
# ╟─c538e51c-8365-4070-b471-b1d4c49a3511
# ╟─2146995a-4493-4b51-8458-a266995a566e
# ╟─cceb8405-6585-4798-b665-a90fb032b1af
# ╟─3d93cc35-b1eb-4efd-b7ba-4279c171692f
# ╟─b7ab44db-4cdb-4b75-a322-c8e997dcbc65
# ╟─c071bc59-7701-4533-9e3f-76f451b2f140
# ╟─ea084619-514a-4fa6-b7ff-ed749b27cccb
# ╟─24156aa7-f981-412e-b05e-52e6eec5c9db
# ╟─f9ca3c6d-59a1-41cc-a955-2a8a39d4a6cc
# ╟─b7503956-1412-4e7d-855b-480d1209414f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
