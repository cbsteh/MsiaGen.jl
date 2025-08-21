
@with_kw mutable struct Wind <: AbstractMetParam
    mean::Vector{Float64} = []
    sd::Vector{Float64} = []
    rlag::Vector{Float64} = []
end


function create_wind(df::AbstractDataFrame)
    colnames = names(df)
    kw = "wind"
    t = Met{Wind}[]
    for g ∈ groupby(df, :year)
        obs = Wind()
        obs.mean = first_row_to_vec(select(g, startswith.(colnames, "mean_$(kw)")))
        obs.sd = first_row_to_vec(select(g, startswith.(colnames, "sd_$(kw)")))
        obs.rlag = first_row_to_vec(select(g, startswith.(colnames, "rlag_$(kw)")))
        push!(t, Met{Wind}(year=g.year[1], obs=obs))
    end
    t
end


function generate!(wind::Met{Wind}; verbose::Bool=true)
    @unpack year, obs, est = wind

    yr_mean = obs.mean[@m 0]
    yr_sd = obs.sd[@m 0]
    yr_rlag = obs.rlag[@m 0]

    thd = [5.0, 5.0, 5.0]
    tgt = [yr_mean, yr_sd, yr_rlag]

    verbose && print_start(year, tgt, thd)

    sz = isleapyear(year) ? 366 : 365
    data = zeros(sz)
    monthdays = days_in_each_month(year)

    generate_monthly_wind!(data, monthdays, obs)

    # CBSTEH:
    data .*= 1.0

    est_mth_mean, est_sd, est_rlag = wind_estimates(data, monthdays)
    est_yr_mean = est_mth_mean[@m 0]

    allok, err, delta = check_errors(thd, tgt, [est_yr_mean, est_sd, est_rlag])

    wind.errors = err
    wind.values = data

    sbm = summarize_by_month
    est.mean = [mean(data), sbm(year, data, mean)...]
    est.sd = [est_sd, sbm(year, data, std)...]
    est.rlag = [est_rlag, sbm(year, data, acf1)...]

    verbose && print_update(allok, wind.errors)
end


# generate daily wind speed for a given month:
function wind_dist!(data, i0, i1, avg, sd, rlag)
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

    sz = i1 - i0
    i2 = i1 - 1
    min_err = 999_999_999.99
    maxrun = 5_000
    nrun = 0
    bOk = false

    while !(min_err <= 2.5) && (nrun < maxrun)
        nrun += 1
        e = rand(wdist, sz) .- avg        # subtract mean from every residual value

        if i0 == 1
            data[i0] = max(0.1, c + rlag * avg + e[i0])
            i0 = 2
        end

        for pos ∈ i0:i2
            data[pos] = max(0.1, c + rlag * data[pos-1] + e[pos-i0+1])
        end

        x = @view data[i0:i2]
        err_avg = 100 * abs((mean(x) - avg) / avg)
        err_sd = 100 * abs((std(x) - sd) / sd)
        err_rlag = 100 * abs((acf1(x) - rlag) / rlag)
        err = max(err_avg, err_sd, err_rlag)

        if err < min_err
            min_err = err
        end
    end
end


function generate_monthly_wind!(data, monthdays, obs)
    t0 = 1
    for i ∈ 1:12
        t1 = t0 + monthdays[i]
        avg = obs.mean[@m i]
        sd = obs.sd[@m i]
        rlag = obs.rlag[@m i]
        wind_dist!(data, t0, t1, avg, sd, rlag)
        t0 = t1
    end
end


function wind_estimates(data, lengths)
    est_mth_mean = mean.(partition_list(data, lengths))
    est_sd = std(data)
    est_rlag = acf1(data)
    est_mth_mean, est_sd, est_rlag
end
