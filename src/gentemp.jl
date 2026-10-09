
@with_kw mutable struct Temp <: AbstractMetParam
    mean::Vector{Float64} = []
    sd::Vector{Float64} = []
    rlag::Vector{Float64} = []
    skew::Vector{Float64} = []
end


create_temp(df::AbstractDataFrame, kw::AbstractString) = create_mets(Temp, df, kw)


"""Skewed normal distribution (for normal/non-extreme skews), written into `e`;
`u` is a work vector of the same length"""
function normal_skew!(e, u, avg, sd, skew)
    sk2_3 = abs(skew) ^ (2 / 3)
    n1 = 0.5 * pi * sk2_3
    n2 = sk2_3 + ((4 - pi) / 2) ^ (2 / 3)
    delta = copysign(sqrt(n1 / n2), skew)     # delta has same sign with skew
    shape = delta / sqrt(1 - delta ^ 2)       # skew <0.995272, so that delta <1

    scale = sqrt(sd ^ 2 / (1 - 2 * delta ^ 2 / pi))
    loc = avg - scale * sqrt(2 / pi) * delta

    randn!(e)
    randn!(u)
    for k ∈ eachindex(e, u)
        z = (u[k] > shape * e[k]) ? -e[k] : e[k]
        e[k] = loc + scale * z
    end
    e
end


"""Extreme skewed normal distribution" (for high skews), written into `e`"""
function high_skew!(e, avg, sd, skew)
    # simulate a skewed normal using an F-distribution:
    #   a) fix the second df2 to 500, then
    #   b) work out the first df1 so that the F-distribution has the desired skew
    sk = abs(skew)  # F-distribution is always skewed right (+ve)
    sk2 = sk^2
    df2 = 500
    d6 = df2 - 6
    d4 = df2 - 4
    d2 = df2 - 2
    a = sqrt(-32 * d4 + sk2 * d6^2)
    b = d2 * (-d6 * sk + a)
    df1 = -b / (2 * a)
    rand!(FDist(df1, df2), e)  # sample the F-distribution

    # transform the sampled F-distribution to have the desired mean and SD
    fmean = mean(e)
    fsd = std(e)
    e .= avg .+ (e .- fmean) .* sd ./ fsd

    if skew < 0
        m = mean(e)
        e .= m .- e  # flip the distribution for negative skew
    end

    # finally, adjust the distribution mean to the specified value
    e .+= avg - mean(e)
end


function skewnorm_rvs!(e, u, avg, sd, skew)
    # max. |skew| shoule be lower at <0.995272 (not 0.99552717) to avoid
    #    calculations later on using imaginary numbers (complex values)
    abs(skew) < 0.995272 ? normal_skew!(e, u, avg, sd, skew) : high_skew!(e, avg, sd, skew)
end


# Thresholds (%) of the whole-year errors shown with `verbose`
fit_thresholds(::Type{Temp}, tgt) = [5.0, 5.0, 10 / max(abs(tgt[3]), 1e-3),
                                     10 / max(abs(tgt[4]), 1e-3)]


# Fit tolerances of a month's generated temperatures: mean (°C), sd
# (fraction of the target), rlag and skew (absolute)
const TEMP_TOL = (mean=0.1, sd=0.025, rlag=0.025, skew=0.05)


# How much warmer (+) or cooler (-) Tmax is on wet days than on dry days of
# the same month (°C): the mean of 23 Malaysian sites (observed daily
# weather, 3-5 years each; site values -0.30 to -0.93, SE 0.04; see
# check_links). Wet days have rain above 0.
const WETDRY_TMAX = -0.6


# Mean of each day of a month, relative to the month's mean, so that wet
# days (`wet`) are `wetdry` warmer than dry days and the month's mean is
# unchanged; and the variance and lag-1 autocovariance these day means add
# to the month. The difference is reduced if it would take more than half
# the month's variance `sd^2`.
function wetdry_shift(wet, wetdry, sd)
    f = mean(wet)                         # share of wet days
    v = f * (1 - f) * wetdry^2            # variance added by the shift
    v > 0.5 * sd^2 && (wetdry *= sqrt(0.5 * sd^2 / v); v = 0.5 * sd^2)
    shift = [w ? (1 - f) * wetdry : -f * wetdry for w ∈ wet]
    shift, v, v * acf1(Float64.(wet))
end


# Generate month i (days `r` of `data`); `prev` is the day before the month.
# With `wet` (each day wet or not), Tmax is `wetdry` (°C) warmer on wet days
# than on dry days; the autoregression then makes up only the rest of the
# month's variance and lag-1 autocovariance.
function generate_month!(obs::Temp, i, data, r, prev; wet=nothing, wetdry=0.0)
    avg = obs.mean[@m i]
    sd = obs.sd[@m i]
    rlag = obs.rlag[@m i]
    skew = obs.skew[@m i]

    shift, v, cv = isnothing(wet) || iszero(wetdry) ? (nothing, 0.0, 0.0) :
                   wetdry_shift(wet, wetdry, sd)
    var_a = sd^2 - v
    rlag_a = clamp((rlag * sd^2 - cv) / var_a, -0.99, 0.99)
    sde = sqrt(var_a * (1 - rlag_a^2))
    c = avg * (1 - rlag_a)
    u = zeros(length(r))

    # each error relative to its tolerance (TEMP_TOL); 1 or less is a fit.
    # Absolute tolerances for mean, rlag and skew, whose targets can be
    # large (mean) or near 0 (rlag, skew); relative for sd.
    score(s) = max(abs(s.mean - avg) / TEMP_TOL.mean,
                   abs(s.sd - sd) / (TEMP_TOL.sd * sd),
                   abs(s.rlag - rlag) / TEMP_TOL.rlag,
                   abs(s.skew - skew) / TEMP_TOL.skew)

    autoregress_month!(data, r, prev, c, rlag_a, e -> skewnorm_rvs!(e, u, 0.0, sde, skew), score;
                       shift=shift)
end
