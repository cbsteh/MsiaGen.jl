
@with_kw mutable struct Temp <: AbstractMetParam
    mean::Vector{Float64} = []
    sd::Vector{Float64} = []
    rlag::Vector{Float64} = []
    skew::Vector{Float64} = []
end


function create_temp(df::AbstractDataFrame, kw::AbstractString)
    colnames = names(df)
    t = Met{Temp}[]
    for g ∈ groupby(df, :year)
        obs = Temp()
        obs.mean = first_row_to_vec(select(g, startswith.(colnames, "mean_$(kw)")))
        obs.sd = first_row_to_vec(select(g, startswith.(colnames, "sd_$(kw)")))
        obs.rlag = first_row_to_vec(select(g, startswith.(colnames, "rlag_$(kw)")))
        obs.skew = first_row_to_vec(select(g, startswith.(colnames, "skew_$(kw)")))
        push!(t, Met{Temp}(year=g.year[1], obs=obs))
    end
    t
end


"""Skewed normal distribution (for normal/non-extreme skews)"""
function normal_skew(avg, sd, skew, sz)
    sk2_3 = abs(skew) ^ (2 / 3)
    n1 = 0.5 * pi * sk2_3
    n2 = sk2_3 + ((4 - pi) / 2) ^ (2 / 3)
    delta = copysign(sqrt(n1 / n2), skew)     # delta has same sign with skew
    shape = delta / sqrt(1 - delta ^ 2)       # skew <0.995272, so that delta <1

    dlt = shape / sqrt(1 + shape ^ 2)
    scale = sqrt(sd ^ 2 / (1 - 2 * dlt ^ 2 / pi))
    loc = avg - scale * sqrt(2 / pi) * dlt

    u1 = rand(Normal(0.0, 1.0), sz)
    u2 = rand(Normal(0.0, 1.0), sz)
    i = findall(u2 .> shape .* u1)
    u1[i] .*= -1.0
    loc .+ scale .* u1
end


"""Extreme skewed normal distribution" (for high skews)"""
function high_skew(avg, sd, skew, sz)
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
    fsample = rand(FDist(df1, df2), sz)  # sample the F-distribution

    # transform the sampled F-distribution to have the desired mean and SD
    fmean = mean(fsample)
    fsd = std(fsample)
    new_dist = avg .+ (fsample .- fmean) .* sd ./ fsd

    if skew < 0
        new_dist = mean(new_dist) .- new_dist  # flip the distribution for negative skew
    end

    # finally, adjust the distribution mean to the specified value
    new_dist .+ (avg .- mean(new_dist))
end


function skewnorm_rvs(avg, sd, skew, sz)
    # max. |skew| shoule be lower at <0.995272 (not 0.99552717) to avoid
    #    calculations later on using imaginary numbers (complex values)
    fn = abs(skew) < 0.995272 ? normal_skew : high_skew
    fn(avg, sd, skew, sz)
end


function generate!(temp::Met{Temp}; verbose::Bool=true)
    @unpack year, obs, est = temp

    yr_mean = obs.mean[@m 0]
    yr_sd = obs.sd[@m 0]
    yr_rlag = obs.rlag[@m 0]
    yr_skew = obs.skew[@m 0]

    thd = [5.0, 5.0, 10 / abs(yr_rlag), 10 / abs(yr_skew)]
    tgt = [yr_mean, yr_sd, yr_rlag, yr_skew]

    verbose && print_start(year, tgt, thd)

    sz = isleapyear(year) ? 366 : 365
    data = zeros(sz)
    monthdays = days_in_each_month(year)

    data = generate_monthly_temp!(data, monthdays, obs, sz, yr_mean)
    est_mth_mean, est_sd, est_rlag, est_skew = temp_estimates(data, monthdays)
    est_yr_mean = est_mth_mean[@m 0]

    est_lst = [est_yr_mean, est_sd, est_rlag, est_skew]
    allok, err, delta = check_errors(thd, tgt, est_lst)

    temp.errors = err
    temp.values = data

    sbm = summarize_by_month
    est.mean = [mean(data), sbm(year, data, mean)...]
    est.sd = [est_sd, sbm(year, data, std)...]
    est.rlag = [est_rlag, sbm(year, data, acf1)...]
    est.skew = [est_skew, sbm(year, data, skewness)...]

    verbose && print_update(allok, temp.errors)
end


function temp_dist!(data, i0, i1, avg, sd, rlag, skew, sz, yr_avg)
    sde = sqrt((sd^2) * (1 - rlag^2))
    c = avg * (1 - rlag)

    sz = i1 - i0
    i2 = i1 - 1
    min_err = 999_999_999.99
    maxrun = 5_000
    nrun = 0
    bOk = false

    while !(min_err <= 2.5) && (nrun < maxrun)
        nrun += 1
        e = skewnorm_rvs(0.0, sde, skew, sz)

        if i0 == 1
            data[i0] = c + rlag * yr_avg + e[i0]
            i0 = 2
        end

        for pos in i0:i2
            data[pos] = c + rlag * data[pos-1] + e[pos-i0+1]
        end

        x = @view data[i0:i2]
        err_avg = 100 * abs((mean(x) - avg) / avg)
        err_sd = 100 * abs((std(x) - sd) / sd)
        err_rlag = 100 * abs((acf1(x) - rlag) / rlag)
        err_skew = 100 * abs((skewness(x) - skew) / skew)
        err = max(err_avg, err_sd, err_rlag, err_skew)

        if err < min_err
            min_err = err
        end
    end
end


function generate_monthly_temp!(data, monthdays, obs, sz, yr_avg)
    t0 = 1
    for i ∈ 1:12
        t1 = t0 + monthdays[i]
        avg = obs.mean[@m i]
        sd = obs.sd[@m i]
        rlag = obs.rlag[@m i]
        skew = obs.skew[@m i]
        temp_dist!(data, t0, t1, avg, sd, rlag, skew, sz, yr_avg)
        t0 = t1
    end
    data
end


function temp_estimates(data, lengths)
    est_mth_mean = mean.(partition_list(data, lengths))
    est_sd = std(data)
    est_rlag = acf1(data)
    est_skew = skewness(data)
    est_mth_mean, est_sd, est_rlag, est_skew
end
