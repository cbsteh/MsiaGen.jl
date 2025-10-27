
@with_kw mutable struct Rain <: AbstractMetParam
    totrain::Vector{Float64} = []
    pww::Vector{Float64} = []
    pwd::Vector{Float64} = []
end


function create_rain(df::AbstractDataFrame)
    colnames = names(df)
    t = Met{Rain}[]
    for g ∈ groupby(df, :year)
        obs = Rain()
        obs.totrain = first_row_to_vec(select(g, startswith.(colnames, "totrain")))
        obs.pww = first_row_to_vec(select(g, startswith.(colnames, "pww")))
        obs.pwd = first_row_to_vec(select(g, startswith.(colnames, "pwd")))
        push!(t, Met{Rain}(year=g.year[1], obs=obs))
    end
    t
end


# return month no. based on day of year
function month(doy::Int, cumdays::AbstractVector)
    i0 = (doy < 31) ? 1 : findlast(cumdays .<= doy)
    i1 = (doy >= 365) ? 12 : findfirst(cumdays .> doy)
    (cumdays[i0] < doy <= cumdays[i1]) ? i1 : i0
end


# Function to count WW and WD days
function collect_rain_counts(amt::AbstractVector)
    nww = nwd = 0
    nw = (amt[1] > 0.0) ? 1 : 0
    N = length(amt)
    for i ∈ 2:N
        w = (amt[i] > 0.0) ? 1 : 0
        nw += w
        ww = (amt[i-1] > 0.0) && (w > 0) ? 1 : 0
        nww += ww
        wd = (amt[i-1] ≈ 0.0) && (w > 0) ? 1 : 0
        nwd += wd
    end
    nd = N - nw
    nw, nd, nww, nwd
end


function prob_WW_WD(amt::AbstractVector)
    nw, nd, nww, nwd = collect_rain_counts(amt)
    pww = (nw > 0) ? nww / nw : 0.0
    pwd = (nd > 0) ? nwd / nd : 0.0
    pw = nw / (nw + nd)
    pww, pwd, pw
end


function rain_stats(amt::AbstractVector)
    pww, pwd, pw = prob_WW_WD(amt)
    n_mths = 12 * length(amt) / 365
    amtx = amt[amt .> 0.0]
    totrain = sum(amtx)
    avg_mth = totrain / n_mths
    (; totrain=totrain, mean=avg_mth, pww=pww, pwd=pwd, pw=pw)
end


function partition_rain(year::Int, dailyrain::AbstractVector)
    monthdays = days_in_each_month(year)
    rain_each_month = partition_list(dailyrain, monthdays)

    totrain = sum.(rain_each_month)     # monthly rainfalls
    pushfirst!(totrain, sum(totrain))   # annual rainfall

    pww_lst = Float64[]
    pwd_lst = Float64[]

    pww, pwd = prob_WW_WD(dailyrain)
    push!(pww_lst, pww)
    push!(pwd_lst, pwd)

    for r ∈ rain_each_month
        pww, pwd = prob_WW_WD(r)
        push!(pww_lst, pww)
        push!(pwd_lst, pwd)
    end

    totrain, pww_lst, pwd_lst
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


function gen_wetdays(sz::Int, totrain, μ)
    x = zeros(sz)
    isapprox(μ, 0) && return x

    min_err = 999_999_999.99
    maxrun = 1_000
    nrun = 0

    while !(min_err <= 2.5) && (nrun < maxrun)
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

    x
end


function distribute_wetdays(sz::Int, x, pww, pwd, pw, rain0)
    min_err = 999_999_999.99
    maxrun = 1_000
    nrun = 0

    szx = length(x)
    finalx = zeros(sz)

    while !(min_err <= 5.0) && (nrun < maxrun)
        nrun += 1
        rs = rand(sz)
        ix = 1
        est_x = zeros(sz)

        w1 = (rain0 < 0.0) ? (rs[1] <= pw) : (rain0 > 0.0)

        for (i, r) ∈ enumerate(rs)
            w0 = w1
            p = w0 ? pww : pwd
            w1 = (r <= p)
            if w1 && (ix <= szx)
                est_x[i] = x[ix]
                ix += 1
            end
        end

        Δx = szx - length(est_x[est_x .> 0.0])
        if Δx > 0
            # not enough wet days; convert some dry days to wet days (at random)
            idx = findall(v->isapprox(v, 0.0), est_x)
            s0 = (length(idx) >= Δx) ? sample(idx, Δx; replace=false) : idx
            balance = length(s0)     # this will be either Δx or length(idx)
            est_x[s0] .= x[end-balance+1:end]
        end

        st = rain_stats(est_x)
        err_pww = 100 * abs(st.pww - pww) / max(0.01, pww)
        err_pwd = 100 * abs(st.pwd - pwd) / max(0.01, pwd)
        err = max(err_pww, err_pwd)

        if err < min_err
            min_err = err
            finalx = est_x[:]
        end
    end

    finalx, min_err
end


function gen_rain_month(sz, totrain, pww, pwd, rain0)
    d = 1 - pww + pwd
    pw = isapprox(d, 0.0) ? 1.0 : pwd / d
    nw = Int(floor(sz * pw))
    μ = (nw > 0) ? totrain / nw : 0.0   # a month may be completely rain-free
    x = gen_wetdays(nw, totrain, μ)
    distribute_wetdays(sz, x, pww, pwd, pw, rain0)
end


function generate!(rain::Met{Rain}; verbose::Bool=true)
    @unpack year, obs, est = rain

    thd = [5.0, 10.0, 10.0]
    tgt = [obs.totrain[@m 0], obs.pww[@m 0], obs.pwd[@m 0]]

    verbose && print_start(year, tgt, thd)

    daysmth = days_in_each_month(year)
    x = [Float64[] for _ ∈ 1:12]

    for i ∈ 1:12
        sz = daysmth[i]
        totrain = obs.totrain[@m i]
        pww = obs.pww[@m i]
        pwd = obs.pwd[@m i]
        rain0 = (i > 1) ? x[i-1][end] : -1.0
        x[i], err = gen_rain_month(sz, totrain, pww, pwd, rain0)
    end

    est_dailyrain = collect(Iterators.flatten(x))

    # determine the fitting errors:
    est_totrain, est_pww, est_pwd = partition_rain(year, est_dailyrain)
    est_val = [est_totrain[@m 0], est_pww[@m 0], est_pwd[@m 0]]

    allok, err, delta = check_errors(thd, tgt, est_val)

    rain.errors = err
    est.totrain = est_totrain
    est.pww = est_pww
    est.pwd = est_pwd
    rain.values = est_dailyrain

    verbose && print_update(allok, rain.errors)
end
