# The weather designer's science: the statistics of every year built from the design
# (start-year values, month changes, end-year targets and the path to them), the
# illustrative days drawn with MsiaGen's own equations, the checks, and the stats file.
# The illustrative days follow MsiaGen's own equations (gentemp.jl, genwind.jl, genrain.jl).
#
# A design is one Dict per part ("tmin", "tmax", "wind", "rain"), each with:
#   "whole"  => the start-year values of the whole year
#   "month"  => the changes of each month (12 values per statistic)
#   "target" => the path ("linear", "accel", "level", "scurve") and the change of each
#               whole-year value by the end year (rain total: in %, others: in its unit)
# The end-year target is the start value plus the change, so changing the start year
# keeps the scenario.

const MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# slider steps from lo to hi (adding 0.0 turns a rounded -0.0 into 0.0)
steps(lo, hi, step) = round.(range(lo, hi; step=step); digits=3) .+ 0.0
snap(r, x) = r[argmin(abs.(r .- x))]

# days in each month of year y, and the month of every day
function calendar(y)
    mdays = [daysinmonth(y, m) for m ∈ 1:12]
    (mdays=mdays, ndays=sum(mdays), month_of=reduce(vcat, [fill(m, mdays[m]) for m ∈ 1:12]))
end

# How a value moves from its start-year value to its end-year target: g(f) is the share
# of the change reached at fraction f of the way from the start year (f = 0) to the end
# year (f = 1)
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

# how far year y is from the start year to the end year (0 in the start year, 1 in the end)
fraction(years) = y -> length(years) == 1 ? 0.0 : (y - first(years)) / (last(years) - first(years))

# ---- reading a daily weather file and its average year ------------------------------------

# Reads a daily weather CSV: an optional latitude line, '#' comment lines, a header, then
# numbers (missing values become NaN).
function parse_weather(text::AbstractString)
    lines = [strip(l) for l ∈ split(text, '\n')]
    lines = [l for l ∈ lines if !isempty(l) && !startswith(l, "#")]
    isempty(lines) && error("the file is empty")
    lat = tryparse(Float64, first(lines))
    isnothing(lat) || (lines = lines[2:end])
    isempty(lines) && error("the file has no header")
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
# pww and pwd over the days that have a next day (all but the last); no dry days (every
# day wet): pwd = 1, as in MsiaGen
function rain_probs(r)
    wet = r .> 0
    prev, next = wet[1:end-1], wet[2:end]
    nw, nd = count(prev), count(.!prev)
    (nw > 0 ? count(prev .& next) / nw : 0.0, nd > 0 ? count(.!prev .& next) / nd : 1.0)
end

# The file's average year: each statistic of each month, from every year, averaged over
# the years. A variable missing from the file is `nothing`.
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

# Fitting the designer's seasonal curves to 12 monthly values, on a 365-day year (so the
# start positions do not depend on the start year)
const CAL365 = calendar(2001)
const MONTH_COS = [mean(cos.(2π .* findall(CAL365.month_of .== m) ./ 365)) for m ∈ 1:12]
const MONTH_SIN = [mean(sin.(2π .* findall(CAL365.month_of .== m) ./ 365)) for m ∈ 1:12]
day_angle(day) = 2π * day / 365
angle_day(θ) = mod(round(Int, θ * 365 / (2π)) - 1, 365) + 1

# monthly means of the curve M + A cos(2π (day - peak) / 365), and the curve that fits 12
# monthly means best
curve_months(M, A, peak) = M .+ A .* (cos(day_angle(peak)) .* MONTH_COS .+
                                      sin(day_angle(peak)) .* MONTH_SIN)
function fit_curve(μ)
    β = [ones(12) MONTH_COS MONTH_SIN] \ μ
    (M=β[1], A=hypot(β[2], β[3]), peak=angle_day(atan(β[3], β[2])))
end
# monthly shares of the rain curve 1 - s cos(2π (day - dry) / 365), and the curve that
# fits 12 monthly totals best
function rain_shares(seas, dry)
    curve = 1 .- (seas / 100) .* cos.(2π .* ((1:365) .- dry) ./ 365)
    [sum(curve[CAL365.month_of .== m]) for m ∈ 1:12] ./ sum(curve)
end
function fit_rain(total)
    β = [ones(12) MONTH_COS MONTH_SIN] \ (total ./ CAL365.mdays)
    a, b = β[2] / β[1], β[3] / β[1]
    (seas=100 * hypot(a, b), dry=angle_day(atan(-b, -a)))
end

# ---- air temperature -------------------------------------------------------------------------

const TVARS = ["tmin", "tmax"]
# Start values of the whole year, used when no weather file is loaded. sd, rlag and skew
# are the medians of all monthly values in the 23 sites' stats files.
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

nomonths(keys) = Dict{String, Any}(k => zeros(12) for k ∈ keys)
target0(keys) = Dict{String, Any}("path" => "linear", (k => 0.0 for k ∈ keys)...)

# The design of tmin or tmax: from the loaded file's average year (whole-year curve fitted
# to its monthly means; each month holds that month's remainder, all snapped to the
# slider steps), or the start values with no month changes. `b`: the file's statistics of
# the variable, or nothing.
function temp_design(v, b)
    keys = ("M", "A", "peak", "sd", "rlag", "skew")
    if isnothing(b)
        d = T_START[v]
        whole = Dict{String, Any}(String(k) => getfield(d, k) for k ∈ propertynames(d))
        month = nomonths(("dm", "dsd", "drl", "dsk"))
    else
        r = T_RANGE
        f = fit_curve(b.μ)
        M, A = snap(r.M, f.M), snap(r.A, f.A)
        sd, rl, sk = snap(r.sd, mean(b.σ)), snap(r.rlag, mean(b.ρ)), snap(r.skew, mean(b.γ))
        whole = Dict{String, Any}("M" => M, "A" => A, "peak" => f.peak, "sd" => sd, "rlag" => rl, "skew" => sk)
        month = Dict{String, Any}("dm" => snap.(Ref(r.dm), b.μ .- curve_months(M, A, f.peak)),
                                  "dsd" => snap.(Ref(r.dsd), b.σ .- sd), "drl" => snap.(Ref(r.drl), b.ρ .- rl),
                                  "dsk" => snap.(Ref(r.dsk), b.γ .- sk))
    end
    Dict{String, Any}("whole" => whole, "month" => month, "target" => target0(keys))
end

num(x) = Float64(x)
months(d, k) = num.(d["month"][k])
# the start value of `k` and its end-year target (start + change)
start_end(d, k) = (w = num(d["whole"][k]); (w, w + num(d["target"][k])))

# Monthly temperature statistics of design `d` (tmin or tmax) in year y, at fraction f of
# the way to the end year. The whole-year values move from the start year towards the
# end-year targets along the chosen path; the year's seasonal curve gives its 12 monthly
# means; the month changes are then added. Values beyond T_LIMITS are capped, and the
# year is flagged.
function t_year_stats(d, y, f)
    g = path_g(d["target"]["path"], f)
    tw(k) = toward(start_end(d, k)..., g)
    M = tw("M")
    A = max(0.0, tw("A"))
    peak = toward_day(start_end(d, "peak")..., g)
    cal = calendar(y)
    curve = M .+ A .* cos.(2π .* ((1:cal.ndays) .- peak) ./ cal.ndays)
    μ = [mean(curve[cal.month_of .== m]) for m ∈ 1:12] .+ months(d, "dm")
    raw = (sd=tw("sd") .+ months(d, "dsd"),
           rlag=tw("rlag") .+ months(d, "drl"),
           skew=tw("skew") .+ months(d, "dsk"))
    capped(k) = round.(clamp.(raw[k], T_LIMITS[k]...); digits=3)
    (year=y, M=M, A=A, peak=peak, μ=μ,
     σ=capped(:sd), ρ=capped(:rlag), γ=capped(:skew),
     mdays=cal.mdays, ndays=cal.ndays,
     clamped=any(any(raw[k] .!= clamp.(raw[k], T_LIMITS[k]...)) for k ∈ keys(raw)))
end

# Skewed residuals with mean 0 and SD `sde`, using MsiaGen's two methods (gentemp.jl:
# normal_skew and high_skew), driven by fixed random numbers.
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

# fixed random numbers, so points move smoothly when a slider changes
t_draws(n) = (rng = MersenneTwister(1234 + n); (u1=randn(rng, 366), u2=randn(rng, 366), uf=rand(rng, 366)))

# One illustrative year of temperature: MsiaGen's lag-1 autoregression, month by month
# (gentemp.jl: temp_dist!), with fixed random numbers.
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

other_tvar(v) = v == "tmin" ? "tmax" : "tmin"

# Skew-normal distribution with a month's mean, sd and skew, as MsiaGen draws daily
# temperatures (|skew| taken as at most 0.99)
function skewnormal(avg, sd, skew)
    sk = clamp(skew, -0.99, 0.99)
    sk2_3 = abs(sk)^(2 / 3)
    δ = copysign(sqrt(0.5π * sk2_3 / (sk2_3 + ((4 - π) / 2)^(2 / 3))), sk)
    ω = sd / sqrt(1 - 2δ^2 / π)
    SkewNormal(avg - ω * sqrt(2 / π) * δ, ω, δ / sqrt(1 - δ^2))
end

# Chance that a day's tmin reaches its tmax, from the month's (mean, sd, skew) of tmax
# `tx` and tmin `tn`, drawn independently as MsiaGen does. Within about 15% of the days
# MsiaGen repairs, at 5 observed sites.
function p_cross(tx, tn)
    dx, dn = skewnormal(tx...), skewnormal(tn...)
    w = 8 * max(tx[2], tn[2])
    t = range(min(tx[1], tn[1]) - w, max(tx[1], tn[1]) + w; length=2001)
    above = reverse(cumsum(reverse(pdf.(dn, t)))) .* step(t)   # chance tmin >= t
    sum(pdf.(dx, t) .* above) * step(t)
end

# Months of all years where tmin may reach tmax on 1% of days or more (year, month,
# expected days), most days first; and the expected days over all years
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
    worst = join(["$(MONTH_NAMES[r.month]) $(r.year) ($(round(r.days; digits=1)) days)"
                  for r ∈ first(rows, 3)], ", ")
    "Tmin may reach tmax on about $(round(total; digits=1)) day(s) over all years, in " *
    "$(length(rows)) month(s) where it is likely on 1% of days or more; most in $(worst). " *
    "MsiaGen swaps tmin and tmax on such days. To avoid it, widen the gap between the tmax " *
    "and tmin means, or lower their sd, in those months."
end

# ---- wind speed ------------------------------------------------------------------------------

# Start values of the whole year, used when no weather file is loaded: medians of the 23
# sites' stats files
const W_START = (M=1.6, A=0.2, peak=45, sd=0.3, rlag=0.2)
# valid limits of each statistic (MsiaGen's wind floor is 0.1 m/s)
const W_LIMITS = (mean=(0.1, 10.0), sd=(0.05, 2.0), rlag=(0.0, 0.95))
const W_FLOOR = 0.1
const W_RANGE = (M=steps(W_LIMITS.mean..., 0.05), A=steps(0.0, 3.0, 0.05), peak=1:365,
                 sd=steps(W_LIMITS.sd..., 0.01), rlag=steps(W_LIMITS.rlag..., 0.01),
                 dm=steps(-2.0, 2.0, 0.05), dsd=steps(-1.0, 1.0, 0.01),
                 drl=steps(-1, 1, 0.1))

function wind_design(b)
    keys = ("M", "A", "peak", "sd", "rlag")
    if isnothing(b)
        d = W_START
        whole = Dict{String, Any}(String(k) => getfield(d, k) for k ∈ propertynames(d))
        month = nomonths(("dm", "dsd", "drl"))
    else
        r = W_RANGE
        f = fit_curve(b.μ)
        M, A = snap(r.M, f.M), snap(r.A, f.A)
        sd, rl = snap(r.sd, mean(b.σ)), snap(r.rlag, mean(b.ρ))
        whole = Dict{String, Any}("M" => M, "A" => A, "peak" => f.peak, "sd" => sd, "rlag" => rl)
        month = Dict{String, Any}("dm" => snap.(Ref(r.dm), b.μ .- curve_months(M, A, f.peak)),
                                  "dsd" => snap.(Ref(r.dsd), b.σ .- sd), "drl" => snap.(Ref(r.drl), b.ρ .- rl))
    end
    Dict{String, Any}("whole" => whole, "month" => month, "target" => target0(keys))
end

# Monthly wind statistics of year y, built like temperature. Values beyond W_LIMITS are
# capped, and the year is flagged.
function w_year_stats(d, y, f)
    g = path_g(d["target"]["path"], f)
    tw(k) = toward(start_end(d, k)..., g)
    M = tw("M")
    A = max(0.0, tw("A"))
    peak = toward_day(start_end(d, "peak")..., g)
    cal = calendar(y)
    curve = M .+ A .* cos.(2π .* ((1:cal.ndays) .- peak) ./ cal.ndays)
    raw = (mean=[mean(curve[cal.month_of .== m]) for m ∈ 1:12] .+ months(d, "dm"),
           sd=tw("sd") .+ months(d, "dsd"),
           rlag=tw("rlag") .+ months(d, "drl"))
    capped(k) = round.(clamp.(raw[k], W_LIMITS[k]...); digits=3)
    (year=y, M=M, A=A, peak=peak,
     μ=capped(:mean), σ=capped(:sd), ρ=capped(:rlag),
     mdays=cal.mdays, ndays=cal.ndays, month_of=cal.month_of,
     clamped=any(any(raw[k] .!= clamp.(raw[k], W_LIMITS[k]...)) for k ∈ keys(raw)))
end

w_draws(n) = rand(MersenneTwister(1234 + n), 366)

# One illustrative year of wind: MsiaGen's wind model, month by month (genwind.jl:
# wind_dist!), with fixed random numbers. Departures come from a Weibull distribution with
# the month's mean and the SD left after the lag-1 term; wind never falls below 0.1 m/s.
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

# Fixed random numbers for the floor check, so it does not flicker
const W_CHECK_U = rand(MersenneTwister(99), 3000)

# Share of days on the 0.1 m/s floor in a month with wind mean μ, sd σ and rlag ρ:
# MsiaGen's recursion run over 3000 days. MsiaGen's search for the best month adds some
# more (up to about 30% more, in tests).
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

# Months of all years where wind may sit on the floor on 1% of days or more
w_floored(w_all) = sort!([(year=s.year, month=m, share=f) for s ∈ w_all for m ∈ 1:12
                          for f ∈ (w_floor_share(s.μ[m], s.σ[m], s.ρ[m]),) if f >= 0.01];
                         by=r -> -r.share)

function w_floor_note(w_all)
    rows = w_floored(w_all)
    isempty(rows) && return nothing
    worst = join(["$(MONTH_NAMES[r.month]) $(r.year) ($(round(Int, 100 * r.share))%)"
                  for r ∈ first(rows, 3)], ", ")
    "In $(length(rows)) month(s), wind may fall to the 0.1 m/s floor on at least 1% of days; " *
    "most in $(worst) (at least that share of days). There, the generated wind sits on the " *
    "floor more often than real wind does, and its mean comes out higher than set. To avoid " *
    "it, lower the sd or raise the mean in those months."
end

# ---- rainfall --------------------------------------------------------------------------------

# Start values of the whole year, used when no weather file is loaded. Annual rainfall,
# pww and pwd are medians of the 23 sites' stats files; seasonality (%) and the driest
# day (day of year) set the seasonal curve.
const R_START = (total=2460.0, seas=30.0, dry=183, pww=0.58, pwd=0.42)
# valid limits: rainfall 0 or more; seasonality 0-90 % (so no month gets negative rain);
# pww and pwd strictly between 0 and 1
const R_LIMITS = (total=(0.0, Inf), seas=(0.0, 90.0), pww=(0.01, 0.99), pwd=(0.01, 0.99))
const R_RANGE = (total=steps(100.0, 6000.0, 10.0), seas=steps(R_LIMITS.seas..., 5.0),
                 dry=1:365, pww=steps(R_LIMITS.pww..., 0.01), pwd=steps(R_LIMITS.pwd..., 0.01),
                 dr=steps(-100.0, 200.0, 1.0), dww=steps(-0.5, 0.5, 0.01),
                 dwd=steps(-0.5, 0.5, 0.01))
# daily rain classes (mm) for the calendar strip
const RAIN_BOUNDS = [0.0, 5.0, 10.0, 20.0, 50.0]

function rain_design(b)
    keys = ("total", "seas", "dry", "pww", "pwd")
    if isnothing(b)
        d = R_START
        whole = Dict{String, Any}(String(k) => getfield(d, k) for k ∈ propertynames(d))
        month = nomonths(("dr", "dww", "dwd"))
    else
        r = R_RANGE
        f = fit_rain(b.total)
        total, seas = snap(r.total, sum(b.total)), snap(r.seas, f.seas)
        pww, pwd = snap(r.pww, mean(b.pww)), snap(r.pwd, mean(b.pwd))
        expected = total .* rain_shares(seas, f.dry)
        whole = Dict{String, Any}("total" => total, "seas" => seas, "dry" => f.dry, "pww" => pww, "pwd" => pwd)
        month = Dict{String, Any}("dr" => snap.(Ref(r.dr), 100 .* (b.total ./ expected .- 1)),
                                  "dww" => snap.(Ref(r.dww), b.pww .- pww), "dwd" => snap.(Ref(r.dwd), b.pwd .- pwd))
    end
    Dict{String, Any}("whole" => whole, "month" => month, "target" => target0(keys))
end

# Expected rain days in a month, as MsiaGen sets them (genrain.jl: gen_rain_month):
# days × pwd / (1 - pww + pwd), rounded, at least 1 if the month has rain and at most the
# days in the month. Average spell lengths of the wet/dry chain.
wet_fraction(pww, pwd) = (d = 1 - pww + pwd; d ≈ 0 ? 1.0 : pwd / d)
rain_days(days, pww, pwd, total) =
    clamp(round(Int, days * wet_fraction(pww, pwd)), total > 0 ? 1 : 0, days)
wet_spell(pww) = 1 / (1 - pww)
dry_spell(pwd) = 1 / pwd

# Monthly rain statistics of year y. The annual rainfall (its end-year target is a change
# in %) is spread over the months by the year's seasonal curve, lowest on the driest day;
# the month changes are then applied. Values beyond R_LIMITS are capped, and the year is
# flagged.
function r_year_stats(d, y, f)
    g = path_g(d["target"]["path"], f)
    w(k) = num(d["whole"][k])
    tw(k) = toward(start_end(d, k)..., g)
    annual = toward(w("total"), w("total") * (1 + num(d["target"]["total"]) / 100), g)
    seas = clamp(tw("seas"), R_LIMITS.seas...)
    dry = toward_day(start_end(d, "dry")..., g)
    cal = calendar(y)
    curve = 1 .- (seas / 100) .* cos.(2π .* ((1:cal.ndays) .- dry) ./ cal.ndays)
    share = [sum(curve[cal.month_of .== m]) for m ∈ 1:12] ./ sum(curve)
    raw = (total=annual .* share .* (1 .+ months(d, "dr") ./ 100),
           pww=tw("pww") .+ months(d, "dww"),
           pwd=tw("pwd") .+ months(d, "dwd"))
    capped(k) = round.(clamp.(raw[k], R_LIMITS[k]...); digits=3)
    total, pww, pwd = capped(:total), capped(:pww), capped(:pwd)
    (year=y, seas=seas, dry=dry, total=total, pww=pww, pwd=pwd, mdays=cal.mdays,
     raindays=rain_days.(cal.mdays, pww, pwd, total),
     wetspell=wet_spell.(pww), dryspell=dry_spell.(pwd),
     clamped=any(any(raw[k] .!= clamp.(raw[k], R_LIMITS[k]...)) for k ∈ keys(raw)))
end

# fixed random numbers, so the strip changes smoothly when a slider moves
function r_draws(n)
    rng = MersenneTwister(1234 + n)
    (chain=rand(rng, 12, 31),     # wet/dry chain
     extra=rand(rng, 12, 31),     # order of dry days turned wet, if needed
     amount=rand(rng, 12, 31),    # rain amounts
     shape=rand(rng, 12))         # Gamma shape of each month
end

# Illustrative daily rain of each month of year s, following MsiaGen (genrain.jl): n =
# expected rain days; n amounts from a Gamma distribution whose shape comes from MsiaGen's
# GEV; wet days placed by the pww/pwd chain (extra chain-wet days stay dry, missing ones
# are taken from dry days). Amounts are scaled so each month's total is exactly the set
# total.
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

# ---- all years, the checks and the stats file -------------------------------------------------

# Statistics of every year of every part
function all_stats(design, years)
    f = fraction(years)
    (t=Dict(v => [t_year_stats(design[v], y, f(y)) for y ∈ years] for v ∈ TVARS),
     w=[w_year_stats(design["wind"], y, f(y)) for y ∈ years],
     r=[r_year_stats(design["rain"], y, f(y)) for y ∈ years])
end

# One MsiaGen stats row per year: tmin, tmax, wind, then rain (the order MsiaGen itself
# writes). Annual (0) columns: temperature and wind mean are day-weighted means of the
# months, rain total is the sum; the other annual values are averages of the 12 months
# (MsiaGen uses them only for its progress messages and the first day of the year).
t_cols(s) = [sum(s.μ .* s.mdays) / s.ndays; s.μ; mean(s.σ); s.σ; mean(s.ρ); s.ρ; mean(s.γ); s.γ]
w_cols(s) = [sum(s.μ .* s.mdays) / s.ndays; s.μ; mean(s.σ); s.σ; mean(s.ρ); s.ρ]
r_cols(s) = [sum(s.total); s.total; mean(s.pww); s.pww; mean(s.pwd); s.pwd]
const STATS_HEADER = ["year";
    ["$(p)_$(v)$(i)" for v ∈ ("tmin", "tmax") for p ∈ ("mean", "sd", "rlag", "skew") for i ∈ 0:12];
    ["$(p)_wind$(i)" for p ∈ ("mean", "sd", "rlag") for i ∈ 0:12];
    ["$(p)$(i)" for p ∈ ("totrain", "pww", "pwd") for i ∈ 0:12]]

# The stats file: latitude line, header, then one row per year
function stats_csv(lat, years, st)
    rows = [[y; t_cols(st.t["tmin"][i]); t_cols(st.t["tmax"][i]); w_cols(st.w[i]); r_cols(st.r[i])]
            for (i, y) ∈ enumerate(years)]
    "$(lat)\n" * join(STATS_HEADER, ",") * "\n" *
        join([join([string(Int(r[1])); string.(round.(r[2:end]; digits=4))], ",") for r ∈ rows], "\n") * "\n"
end

# Years where values were capped (per part)
capped_years(st) = (temperature=sort(unique([s.year for v ∈ TVARS for s ∈ st.t[v] if s.clamped])),
                    wind=[s.year for s ∈ st.w if s.clamped],
                    rain=[s.year for s ∈ st.r if s.clamped])

# The checks shown before the stats file: capped values, tmin reaching tmax, wind on its floor
function value_checks(st)
    notes = String[]
    for (part, ys) ∈ pairs(capped_years(st))
        isempty(ys) || push!(notes, "Some $(part) values went past their limits in $(length(ys)) " *
                                    "year(s), $(first(ys))–$(last(ys)). They were capped at the limit.")
    end
    for n ∈ (t_cross_note(st.t), w_floor_note(st.w))
        isnothing(n) || push!(notes, n)
    end
    notes
end
