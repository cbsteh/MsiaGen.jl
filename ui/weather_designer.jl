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
    using CairoMakie, Dates, Distributions, HypertextLiteral, PlutoUI
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
    const BLUE = colorant"#2a78d6"
    const DARK_BLUE = colorant"#104281"
    const INK2 = colorant"#52514e"
    const GRID = colorant"#e1e0d9"
    const MUTED = colorant"#898781"
    const SURFACE = colorant"#fcfcfb"
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

    # Axis styling shared by all charts
    chart_style = (titlealign=:left, topspinevisible=false, rightspinevisible=false,
               leftspinecolor=GRID, bottomspinecolor=GRID,
               xticklabelcolor=INK2, yticklabelcolor=INK2,
               xlabelcolor=INK2, ylabelcolor=INK2)
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
# Air temperature

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
@htl("""<p><b>Select air temperature type:</b> <span id="tvar-select">$(@bind tvar Select(["tmin" => "Tmin", "tmax" => "Tmax"]))</span></p>""")

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
Markdown.parse("## 1. ($(tvar)) Set start year ($(start_year))")

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
## 2. ($(tvar)) Each month
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
    ## 3. ($(tvar)) Set last year ($(e))

    Set the targets for $(e), and the path the values take from $(start_year) to get there; the years between follow that path, in all months. Targets start at the $(start_year) values, so changing section 1 resets them. The preview below shows $(e).
    """)
end

# ╔═╡ feb14a42-bf45-47f8-a172-bfad503f1db9
# End-year targets of the whole-year values, and the path to them. Targets
# start at the start-year values (section 1).
@bind ttarget PlutoUI.combine() do Child
    parts = map(enumerate(TVARS)) do (i, v)
        r = T_RANGE
        now(k) = twhole[Symbol("$(v)_$(k)")]
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
    function t_year_stats(v, y)
        w(k) = twhole[Symbol("$(v)_$(k)")]
        t(k) = ttarget[Symbol("$(v)_$(k)")]
        mc(k) = [tmonthly[Symbol("$(v)_$(k)$(m)")] for m ∈ 1:12]
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
    function t_chart(s, pts, v, title; ref=nothing)
        mdays, ndays = s.mdays, s.ndays
        starts = cumsum([1; mdays[1:end-1]])
        ends = starts .+ mdays .- 1
        fig = Figure(size=(1000, 460), backgroundcolor=SURFACE, fontsize=16)
        ax = Axis(fig[1, 1]; ylabel="$(v) (°C)", xticks=(starts .+ mdays ./ 2, MONTHS),
                  title=title, xgridvisible=false, ygridcolor=GRID, chart_style...)
        vlines!(ax, starts[2:end] .- 0.5; color=GRID)
        for m ∈ 1:12
            xs = [starts[m] - 0.5, ends[m] + 0.5]
            band!(ax, xs, fill(s.μ[m] - s.σ[m], 2), fill(s.μ[m] + s.σ[m], 2); color=(BLUE, 0.12))
            lines!(ax, xs, fill(s.μ[m], 2); color=DARK_BLUE, linewidth=2.5)
            isnothing(ref) || lines!(ax, xs, fill(ref[m], 2); color=MUTED, linestyle=:dash,
                                     linewidth=2)
        end
        lines!(ax, 1:ndays, pts; color=(BLUE, 0.35), linewidth=0.8)
        scatter!(ax, 1:ndays, pts; color=BLUE, markersize=5)
        xlims!(ax, 0.5, ndays + 0.5)
        elems = Any[MarkerElement(color=BLUE, marker=:circle, markersize=8),
                    LineElement(color=DARK_BLUE, linewidth=2.5), PolyElement(color=(BLUE, 0.12))]
        labels = ["illustrative day", "monthly mean", "±1 sd"]
        isnothing(ref) || (push!(elems, LineElement(color=MUTED, linestyle=:dash, linewidth=2));
                           push!(labels, "start-year monthly mean"))
        axislegend(ax, elems, labels; position=:rb, orientation=:horizontal,
                   framevisible=false, labelcolor=INK2)
        fig
    end
end;

# ╔═╡ 8fc48e5a-3ca5-4279-bf20-b06a237784fa
# temperature statistics of every year, for tmax and tmin
t_all = Dict(v => [t_year_stats(v, y) for y ∈ years] for v ∈ TVARS);

# ╔═╡ edfd7811-cdc2-4a9d-bdc3-99fad649447c
# the start year of the chosen variable, shown in the chart and the table
begin
    t_shown = first(t_all[tvar])
    t_points = t_sample_year(t_shown, t_draws)
end;

# ╔═╡ a88f601f-0b2a-4cad-8818-2054f575d03c
t_chart(t_shown, t_points, tvar,
        "$(tvar) for $(t_shown.year): illustrative days around your monthly statistics") |> WideCell

# ╔═╡ d591c92f-0560-410d-b8d7-f431b6e45396
let s = t_shown, v = tvar
    cell(value, key, m) = "$(value) ($(fmt_change(tmonthly[Symbol("$(v)_$(key)$(m)")])))"
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
    """)
end

# ╔═╡ 86c421d6-d3a0-421f-a6cb-8a856411a21a
# Preview of the end year, with the start year's monthly means for reference
let e = last(t_all[tvar])
    t_chart(e, t_sample_year(e, t_draws), tvar,
            "Preview: $(tvar) for $(e.year), the end year"; ref=t_shown.μ) |> WideCell
end

# ╔═╡ 7b58bb8a-4f4d-49e1-90dc-04b31f3529da
md"""
# Wind speed

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
Markdown.parse("## 1. Set start year ($(start_year))")

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
## 2. Each month
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
    ## 3. Set last year ($(e))

    Set the targets for $(e), and the path the values take from $(start_year) to get there; the years between follow that path, in all months. Targets start at the $(start_year) values, so changing section 1 resets them. The preview below shows $(e).
    """)
end

# ╔═╡ 1f8a6cf2-8e88-4ec3-89e0-357e8e76112f
# End-year targets of the whole-year values, and the path to them. Targets
# start at the start-year values (section 1).
@bind wtarget PlutoUI.combine() do Child
    r = W_RANGE
    now(k) = wwhole[k]
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

# ╔═╡ 561fa940-4e23-4362-9600-ff88dada2f33
begin
    # Monthly wind statistics of year y, built like temperature: whole-year
    # values moved towards the end-year targets along the chosen path, a
    # seasonal curve for the monthly means, then the month sliders. Values
    # beyond W_LIMITS are capped, and the year is flagged.
    function w_year_stats(y)
        w, t = wwhole, wtarget
        mc(k) = [wmonthly[Symbol("$(k)$(m)")] for m ∈ 1:12]
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
    function w_chart(s, pts, title; ref=nothing)
        mdays, ndays = s.mdays, s.ndays
        starts = cumsum([1; mdays[1:end-1]])
        ends = starts .+ mdays .- 1
        fig = Figure(size=(1000, 500), backgroundcolor=SURFACE, fontsize=16)
        ax = Axis(fig[1, 1]; ylabel="wind speed (m/s)", xticks=(starts .+ mdays ./ 2, MONTHS),
                  title=title, xgridvisible=false, ygridcolor=GRID, chart_style...)
        vlines!(ax, starts[2:end] .- 0.5; color=GRID)
        for m ∈ 1:12
            xs = [starts[m] - 0.5, ends[m] + 0.5]
            band!(ax, xs, fill(max(0.0, s.μ[m] - s.σ[m]), 2), fill(s.μ[m] + s.σ[m], 2);
                  color=(BLUE, 0.12))
            lines!(ax, xs, fill(s.μ[m], 2); color=DARK_BLUE, linewidth=2.5)
            isnothing(ref) || lines!(ax, xs, fill(ref[m], 2); color=MUTED, linestyle=:dash,
                                     linewidth=2)
        end
        hlines!(ax, [W_FLOOR]; color=MUTED, linestyle=:dot, linewidth=1)
        lines!(ax, 1:ndays, pts; color=(BLUE, 0.35), linewidth=0.8)
        scatter!(ax, 1:ndays, pts; color=BLUE, markersize=5)
        xlims!(ax, 0.5, ndays + 0.5)
        ylims!(ax, 0, nothing)
        elems = Any[MarkerElement(color=BLUE, marker=:circle, markersize=8),
                    LineElement(color=DARK_BLUE, linewidth=2.5), PolyElement(color=(BLUE, 0.12)),
                    LineElement(color=MUTED, linestyle=:dot)]
        labels = ["illustrative day", "monthly mean", "±1 sd", "0.1 m/s floor"]
        isnothing(ref) || (push!(elems, LineElement(color=MUTED, linestyle=:dash, linewidth=2));
                           push!(labels, "start-year monthly mean"))
        Legend(fig[2, 1], elems, labels; orientation=:horizontal, framevisible=false,
               labelcolor=INK2, halign=:right)
        fig
    end
end;

# ╔═╡ 1ed9263f-1d99-4105-beb2-5a5a6eaa477f
w_all = [w_year_stats(y) for y ∈ years];

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
let s = w_shown
    cell(value, key, m) = "$(value) ($(fmt_change(wmonthly[Symbol("$(key)$(m)")])))"
    rows = map(1:12) do m
        @htl("""<tr><td><b>$(MONTHS[m])</b></td>
                <td>$(cell(round(s.μ[m]; digits=2), "dm", m))</td>
                <td>$(cell(s.σ[m], "sd", m))</td>
                <td>$(cell(s.ρ[m], "rl", m))</td></tr>""")
    end
    floored = [m for m ∈ 1:12 if any(w_points[s.month_of .== m] .<= W_FLOOR)]
    note = isempty(floored) ? "" :
        @htl("""<p style="color: #b52f2f;"><b>Note:</b> some illustrative days sit on the
             0.1 m/s floor in $(join(MONTHS[floored], ", ")). There, the generated mean will
             come out higher than set; lower the sd or raise the mean.</p>""")
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

# ╔═╡ 4f563fdd-9a55-4f0c-a86e-90135154a104
# Preview of the end year, with the start year's monthly means for reference
let e = last(w_all)
    w_chart(e, w_sample_year(e, w_draws), "Preview: wind speed for $(e.year), the end year";
            ref=w_shown.μ) |> WideCell
end

# ╔═╡ 6bbc7b6c-e62d-4244-82e1-422cb6c002c3
md"""
# Rainfall

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
    const RAIN_COLORS = [colorant"#f0efec", colorant"#cde2fb", colorant"#86b6ef",
                         colorant"#3987e5", colorant"#1c5cab", colorant"#0d366b"]
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
Markdown.parse("## 1. Set start year ($(start_year))")

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
## 2. Each month
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
    ## 3. Set last year ($(e))

    Set the targets for $(e), and the path the values take from $(start_year) to get there; the years between follow that path, in all months. Targets start at the $(start_year) values, so changing section 1 resets them. The preview below shows $(e).
    """)
end

# ╔═╡ ad4f0f49-552f-4e2a-9ccc-2ae81c2e6584
# End-year targets of the whole-year values, and the path to them. Targets
# start at the start-year values (section 1).
@bind rtarget PlutoUI.combine() do Child
    r = R_RANGE
    now(k) = rwhole[k]
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
    function r_year_stats(y)
        w, t = rwhole, rtarget
        mc(k) = [rmonthly[Symbol("$(k)$(m)")] for m ∈ 1:12]
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
    # bars as a reference
    function r_chart(s, strip, title; ref=nothing)
        z = fill(NaN, 31, 12)
        for m ∈ 1:12, d ∈ 1:s.mdays[m]
            r = strip[m][d]
            z[d, m] = r > 0 ? searchsortedlast(RAIN_BOUNDS, r) : 0
        end
        fig = Figure(size=(1100, 480), backgroundcolor=SURFACE, fontsize=16)
        ax1 = Axis(fig[1, 1]; yreversed=true, yticks=(1:12, MONTHS), xticks=[1; 5:5:30],
                   xlabel="day of month", xgridvisible=false, ygridvisible=false,
                   title=title, chart_style...)
        heatmap!(ax1, 1:31, 1:12, z; colormap=cgrad(RAIN_COLORS, 6; categorical=true),
                 colorrange=(-0.5, 5.5), nan_color=:transparent)
        vlines!(ax1, 1.5:1:30.5; color=SURFACE, linewidth=2)
        hlines!(ax1, 1.5:1:11.5; color=SURFACE, linewidth=2)
        ax2 = Axis(fig[1, 2]; yreversed=true, yticks=(1:12, MONTHS), xlabel="rainfall (mm)",
                   ygridvisible=false, xgridcolor=GRID, yticklabelsvisible=false,
                   title="Monthly rainfall; label: mm, rain days", chart_style...)
        barplot!(ax2, 1:12, s.total; direction=:x, color=BLUE, gap=0.3)
        text!(ax2, s.total, 1:12;
              text=["$(round(Int, s.total[m]))  $(s.raindays[m]) d" for m ∈ 1:12],
              align=(:left, :center), offset=(6, 0), color=INK2, fontsize=14)
        top = maximum(s.total)
        if !isnothing(ref)
            # a short vertical mark at each month's start-year total
            linesegments!(ax2, [Point2f(ref[m], m - 0.4) => Point2f(ref[m], m + 0.4) for m ∈ 1:12];
                          color=INK2, linewidth=2.5)
            top = max(top, maximum(ref))
        end
        xlims!(ax2, 0, 1.6 * max(1.0, top))
        linkyaxes!(ax1, ax2)
        ylims!(ax1, 12.5, 0.5)
        colsize!(fig.layout, 1, Relative(0.68))
        colgap!(fig.layout, 12)
        # legends under the strip only, so the bar column keeps its width
        leg = GridLayout(fig[2, 1]; halign=:left)
        Legend(leg[1, 1], [PolyElement(color=c, strokecolor=GRID, strokewidth=0.5) for c ∈ RAIN_COLORS],
               ["dry", "< 5", "5–10", "10–20", "20–50", "≥ 50"], "mm/day";
               orientation=:horizontal, titleposition=:left, framevisible=false,
               labelcolor=INK2, titlecolor=INK2)
        isnothing(ref) || Legend(leg[1, 2], [LineElement(color=INK2, linewidth=2.5)],
                                 ["start-year total"]; framevisible=false, labelcolor=INK2)
        fig
    end
end;

# ╔═╡ 2146995a-4493-4b51-8458-a266995a566e
begin
    r_all = [r_year_stats(y) for y ∈ years]
    r_shown = first(r_all)     # the start year, shown in the chart and table
end;

# ╔═╡ a4d95768-97e4-4f28-9244-26918a4ba23f
let s = r_shown
    ch(key, m) = rmonthly[Symbol("$(key)$(m)")]
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
# Output for MsiaGen

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
    stats_rows = [[y; t_cols(t_all["tmin"][i]); t_cols(t_all["tmax"][i]);
                   w_cols(w_all[i]); r_cols(r_all[i])] for (i, y) ∈ enumerate(years)]
end;

# ╔═╡ ea084619-514a-4fa6-b7ff-ed749b27cccb
csv_text = "$(lat)\n" * join(stats_header, ",") * "\n" *
           join([join([string(Int(r[1])); string.(round.(r[2:end]; digits=4))], ",")
                 for r ∈ stats_rows], "\n") * "\n";

# ╔═╡ 24156aa7-f981-412e-b05e-52e6eec5c9db
# download file name: <site>-stats.csv, the name MsiaGen reads (spaces in
# the site name become hyphens)
csv_name = "$(replace(strip(site), r"\s+" => "-"))-stats.csv";

# ╔═╡ f9ca3c6d-59a1-41cc-a955-2a8a39d4a6cc
# Checks: years where values were capped, and months where the tmin mean
# reaches the tmax mean
begin
    capped_years = (temperature=sort(unique([s.year for v ∈ TVARS for s ∈ t_all[v] if s.clamped])),
                    wind=[s.year for s ∈ w_all if s.clamped],
                    rain=[s.year for s ∈ r_all if s.clamped])
    tmin_over_tmax = [(y, MONTHS[m]) for (i, y) ∈ enumerate(years) for m ∈ 1:12
                      if t_all["tmin"][i].μ[m] >= t_all["tmax"][i].μ[m]]
end;

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
    notes = []
    for (part, ys) ∈ pairs(capped_years)
        isempty(ys) || push!(notes, @htl("""<p style="color: #b52f2f;"><b>Note:</b> some
            $(part) values went past their limits in $(length(ys)) year(s), $(first(ys))–$(last(ys)).
            They were capped at the limit.</p>"""))
    end
    isempty(tmin_over_tmax) || push!(notes, @htl("""<p style="color: #b52f2f;"><b>Warning:</b>
        the tmin mean reaches the tmax mean in $(length(tmin_over_tmax)) month(s), first in
        $(tmin_over_tmax[1][2]) $(tmin_over_tmax[1][1]). Raise tmax or lower tmin there.</p>"""))
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CairoMakie = "~0.15.15"
Distributions = "~0.25.131"
HypertextLiteral = "~0.9.5"
PlutoUI = "~0.7.73"
SpecialFunctions = "~2.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.0"
manifest_format = "2.1"
project_hash = "996b144859bdccb69fc33755b64e6be8d5e817f1"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

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

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7c2c19b5a26e601634bf718490b89d59685f122e"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.1"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
registries = "General"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
registries = "General"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
registries = "General"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
registries = "General"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
registries = "General"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "1cda0b7d5abfc95357dae18aca934d401f7869ad"
registries = "General"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "7b841680738c19948120f6e4cf8d0200518bd564"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.8+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
registries = "General"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
registries = "General"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "61761f58648aa7217445f24f841839b78c712232"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.3"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "291665b547f137df070e4dd83e432b5fee8cc4a0"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.2"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

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

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "7bc84b769c1d384315e7b5c4ac03a6c303e6cf35"
registries = "General"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.8"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
registries = "General"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

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

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "4ac548adcad90c1d5d677af13568a748af4c952b"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
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

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
registries = "General"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
registries = "General"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2bfb1e047e2ad0a5ca94365340bde8005d637568"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.4+0"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d9d3cd382f2c999f684e6bfa4039fcc16c153b6f"
registries = "General"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "9.0.2+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
registries = "General"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
registries = "General"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
registries = "General"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
registries = "General"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "086b5fbd032baf544678cc15b52b015e4c4aceb8"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.1"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
registries = "General"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
registries = "General"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
registries = "General"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "ec46c5825710fa1a15d468acb2d93cc939a7a5fe"
registries = "General"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.13"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
registries = "General"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3efbbc027441271444dcd0c0a46f2d119dc4329"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "6.1.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "090526e65de8f69648ac156daae153de8b56df62"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.88.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "ef70da5e123a06a29e2d6ddff0f09985bc226491"
registries = "General"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "9d9531a9cb63a9edc33836414e82a07e81710de2"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.14004.0+0"

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

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "f0f005f997dfb8c5fe23920d99458a9619873893"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.10"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "1c531bf0f8a5c60a340926e058fd3f209b5eef5d"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.12"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticMakieExt = "Makie"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "0db6aea5b64caa1c47e9fdd27394f0296f05e8bf"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.15"

    [deps.IntervalSets.extensions]
    IntervalSetsMakieExt = "Makie"
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

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

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
registries = "General"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
registries = "General"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
registries = "General"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39bca05343661c347aae0bca57a5994a0bf4f08d"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.2.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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
version = "1.11.103+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

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

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "5f6f5d1b1fb7ff98c9a083bbfd4c661a9808e758"
registries = "General"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.15"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
registries = "General"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
registries = "General"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
registries = "General"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e9ae72527609bb66882d38ad009ffd47ca5713fb"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.16+0"

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

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
registries = "General"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

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

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
registries = "General"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
registries = "General"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "ba0dc8a8a67cacac4842631f960c046e4e563675"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.8"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "ddcc36fd83fb4cd459a53c891981697bde0d3de6"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
registries = "General"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

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

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
registries = "General"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

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

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "57aa595158717ef165e6f5ab639fe2e3178c0a2b"
registries = "General"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
registries = "General"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2c79185e5261b159474903598c571ddd9f12ff7b"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
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
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "39e70e0ab5d7f89833a62ab7c79df15d4fc417c1"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.22"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "ad8002667372439f2e3611cfd14097e03fa4bccd"
registries = "General"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.3"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
registries = "General"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
registries = "General"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
registries = "General"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
registries = "General"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
registries = "General"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1210ba774d3427387d307bf1f416d699b7c39417"
registries = "General"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.15.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cb007192783c56d8249db4cf0e3495001edfe414"
registries = "General"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.5+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "28e57478e8a160d346a19c28b3fffb9273bcc9c2"
registries = "General"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.134+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
registries = "General"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "32781be40fe86735af02eae0dd22754e4b5f779d"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "e067c8bae65bb40866552296a2d5b4dc65b8928f"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.80.702+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
registries = "General"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
registries = "General"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "52d3b9475133c3bc8c0a7f90f18bfc8cdb5a443a"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
registries = "General"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
registries = "General"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

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
