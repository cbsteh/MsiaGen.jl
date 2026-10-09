
@with_kw mutable struct Rain <: AbstractMetParam
    totrain::Vector{Float64} = []
    pww::Vector{Float64} = []
    pwd::Vector{Float64} = []
end


create_rain(df::AbstractDataFrame) = create_mets(Rain, df)


# Counts of the days that have a next day (all but the last): wet days
# (nw), dry days (nd), and wet days after a wet day (nww) and a dry day (nwd)
function collect_rain_counts(amt::AbstractVector)
    nw = nd = nww = nwd = 0
    for i ∈ firstindex(amt):lastindex(amt)-1
        wet_next = amt[i+1] > 0.0
        if amt[i] > 0.0
            nw += 1
            nww += wet_next
        else
            nd += 1
            nwd += wet_next
        end
    end
    nw, nd, nww, nwd
end


# Chance of a wet day after a wet day (pww) and after a dry day (pwd). With
# no wet days pww is 0; with no dry days (every day wet) pwd is 1, so the
# month's share of wet days, pwd / (1 - pww + pwd), stays close to 1.
function prob_WW_WD(amt::AbstractVector)
    nw, nd, nww, nwd = collect_rain_counts(amt)
    pww = (nw > 0) ? nww / nw : 0.0
    pwd = (nd > 0) ? nwd / nd : 1.0
    pww, pwd
end


function partition_rain(year::Int, dailyrain::AbstractVector)
    parts = year_and_months(dailyrain, year)     # whole year, then each month
    probs = prob_WW_WD.(parts)
    totrain = sum.(parts[2:end])        # monthly rainfalls
    pushfirst!(totrain, sum(totrain))   # annual rainfall
    totrain, first.(probs), last.(probs)
end


function rand_θ(μ)
    # loc, scale, shape = 0.53231, 0.17817, 0.14581
    loc, scale, shape = 0.50, 0.17, 0.14
    gev = GeneralizedExtremeValue(loc, scale, shape)
    k = -99.9
    while k <= 0
        k = quantile(gev, rand())
    end
    θ = μ / k
    k, θ
end


# Fit tolerances (%) of a month's generated rain: the total, and pww and pwd
const RAIN_TOL = (totrain=2.5, pw=5.0)


# Rain of the month's `sz` wet days, and its fit score (error of the total
# relative to its tolerance; 1 or less is a fit)
function gen_wetdays(sz::Int, totrain, μ)
    x = zeros(sz)
    isapprox(μ, 0) && return x, 0.0

    min_err = 999_999_999.99
    maxrun = 1_000
    nrun = 0

    while !(min_err <= RAIN_TOL.totrain) && (nrun < maxrun)
        nrun += 1
        k, θ = rand_θ(μ)
        est_x = quantile.(Gamma(k, θ), rand(sz))
        est_totrain = sum(est_x)
        err = 100 * abs(est_totrain - totrain) / max(0.01, totrain)

        if err < min_err
            min_err = err
            x = est_x
        end
    end

    x, min_err / RAIN_TOL.totrain
end


# The month's daily rain, the wet days `x` placed by the wet/dry chain, and
# its fit score (largest error of pww and pwd relative to its tolerance)
function distribute_wetdays(sz::Int, x, pww, pwd, pw, rain0)
    min_err = 999_999_999.99
    maxrun = 1_000
    nrun = 0

    szx = length(x)
    finalx = zeros(sz)
    rs = zeros(sz)
    est_x = zeros(sz)

    while !(min_err <= RAIN_TOL.pw) && (nrun < maxrun)
        nrun += 1
        rand!(rs)
        fill!(est_x, 0.0)
        ix = 1

        # wet or dry on the day before the month: as given by `rain0`, or
        # with no previous day, wet with the month's chance of a wet day
        w1 = (rain0 < 0.0) ? (rand() <= pw) : (rain0 > 0.0)

        for (i, r) ∈ enumerate(rs)
            w0 = w1
            p = w0 ? pww : pwd
            w1 = (r <= p)
            if w1 && (ix <= szx)
                est_x[i] = x[ix]
                ix += 1
            end
        end

        Δx = szx - count(>(0.0), est_x)
        if Δx > 0
            # not enough wet days; convert some dry days to wet days (at random)
            idx = findall(v->isapprox(v, 0.0), est_x)
            s0 = (length(idx) >= Δx) ? sample(idx, Δx; replace=false) : idx
            balance = length(s0)     # this will be either Δx or length(idx)
            est_x[s0] .= x[end-balance+1:end]
        end

        est_pww, est_pwd = prob_WW_WD(est_x)
        err_pww = 100 * abs(est_pww - pww) / max(0.01, pww)
        err_pwd = 100 * abs(est_pwd - pwd) / max(0.01, pwd)
        err = max(err_pww, err_pwd)

        if err < min_err
            min_err = err
            copyto!(finalx, est_x)
        end
    end

    finalx, min_err / RAIN_TOL.pw
end


function gen_rain_month(sz, totrain, pww, pwd, rain0)
    d = 1 - pww + pwd
    pw = isapprox(d, 0.0) ? 1.0 : pwd / d     # long-run share of wet days
    # wet days: the expected number, rounded; at least 1 if the month has
    # rain, and never more than the days in the month
    nw = clamp(round(Int, sz * pw), (totrain > 0.0) ? 1 : 0, sz)
    μ = (nw > 0) ? totrain / nw : 0.0   # a month may be completely rain-free
    x, score_total = gen_wetdays(nw, totrain, μ)
    rain, score_pw = distribute_wetdays(sz, x, pww, pwd, pw, rain0)
    rain, max(score_total, score_pw)
end


# `prev`: rain on the day before 1 January (the previous year's 31
# December), or `nothing` to start from the long-run chance of a wet day
function generate!(rain::Met{Rain}; verbose::Bool=true, prev=nothing)
    @unpack year, obs = rain

    thd = [5.0, 10.0, 10.0]
    tgt = [obs.totrain[@m 0], obs.pww[@m 0], obs.pwd[@m 0]]

    verbose && print_start(year, tgt, thd)

    daysmth = days_in_each_month(year)
    x = [Float64[] for _ ∈ 1:12]
    scores = zeros(12)

    for i ∈ 1:12
        sz = daysmth[i]
        totrain = obs.totrain[@m i]
        pww = obs.pww[@m i]
        pwd = obs.pwd[@m i]
        rain0 = (i > 1) ? x[i-1][end] : (isnothing(prev) ? -1.0 : prev)
        x[i], scores[i] = gen_rain_month(sz, totrain, pww, pwd, rain0)
    end

    est_dailyrain = reduce(vcat, x)

    # determine the fitting errors (whole year):
    est_totrain, est_pww, est_pwd = partition_rain(year, est_dailyrain)
    allok, err = check_errors(thd, tgt, [est_totrain[@m 0], est_pww[@m 0], est_pwd[@m 0]])

    rain.errors = err
    rain.values = est_dailyrain
    rain.scores = scores

    verbose && print_update(allok, rain.errors)
end
