
@with_kw mutable struct Wind <: AbstractMetParam
    mean::Vector{Float64} = []
    sd::Vector{Float64} = []
    rlag::Vector{Float64} = []
end


create_wind(df::AbstractDataFrame) = create_mets(Wind, df, "wind")


# Thresholds (%) of the whole-year errors shown with `verbose`
fit_thresholds(::Type{Wind}, _) =[5.0, 5.0, 5.0]


# Fit tolerances of a month's generated wind speeds: mean (m/s), sd
# (fraction of the target), rlag (absolute)
const WIND_TOL = (mean=0.02, sd=0.025, rlag=0.025)


# Generate month i (days `r` of `data`); `prev` is the day before the month
function generate_month!(obs::Wind, i, data, r, prev)
    avg = obs.mean[@m i]
    sd = obs.sd[@m i]
    rlag = obs.rlag[@m i]

    # 1. create a Weibull distribution for the autoregression residuals, BUT
    #   the Weibull distribution is only for +ve values. To have -ve residuals:
    #   a) set the mean of the Weibull distribution equal to the mean data,
    #   b) set the SD of the Weibull distribution equal to the corrected SD of data,
    #   c) determine the shape and scale Weibull parameters, and
    #   d) finally, subtract every Weibull residual by the mean data, so the mean
    #      of all residuals is zero. The Weibull distribution of residuals now have
    #      both +ve and -ve values.
    sde = sqrt((sd^2) * (1 - rlag^2))     # SD corrected for autoregression lag 1
    shape = (sde / avg) ^ -1.086          # shape parameter
    scale = avg / gamma(1 + 1 / shape)    # scale parameter
    wdist = Weibull(shape, scale)         # Weibull residuals
    # 2. autoregression lag 1 equation:
    c = avg * (1 - rlag)  # constant for the autoregression lag 1 equation

    draw!(e) = (rand!(wdist, e); e .-= avg)   # subtract mean from every residual value

    # each error relative to its tolerance (WIND_TOL); 1 or less is a fit.
    # Absolute tolerances for mean and rlag (rlag can be near 0);
    # relative for sd.
    score(s) = max(abs(s.mean - avg) / WIND_TOL.mean,
                   abs(s.sd - sd) / (WIND_TOL.sd * sd),
                   abs(s.rlag - rlag) / WIND_TOL.rlag)

    # wind never falls below 0.1 m/s
    autoregress_month!(data, r, prev, c, rlag, draw!, score; lo=0.1)
end
