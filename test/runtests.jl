using MsiaGen
using CSV
using DataFrames
using Dates
using Random
using Test

const M = MsiaGen


# Made-up daily weather of whole years (the same series as the precompile run)
function test_obs(years)
    dates = Date(first(years)):Day(1):Date(last(years), 12, 31)
    d = 1:length(dates)
    DataFrame(year=Dates.year.(dates), month=Dates.month.(dates), day=Dates.day.(dates),
              tmin=23.0 .+ sin.(d ./ 7) .+ 0.3 .* cos.(d .* 1.3),
              tmax=32.0 .+ sin.(d ./ 5) .+ 0.5 .* cos.(d .* 1.7),
              wind=1.5 .+ 0.3 .* sin.(d ./ 3) .+ 0.2 .* cos.(d .* 2.1),
              rain=[i % 3 == 0 ? 2.0 + i % 11 : 0.0 for i ∈ d])
end

generate(stats; seed=1) = (Random.seed!(seed); M.generate_mets(stats; verbose=false))


@testset "MsiaGen" begin

@testset "statistics of observed weather: dates" begin
    obs = test_obs(2003:2004)
    stats = M.weather_stats(obs)
    @test stats.year == [2003, 2004]

    # rows in any order give the same statistics
    shuffled = obs[randperm(MersenneTwister(1), nrow(obs)), :]
    @test isequal(M.weather_stats(shuffled), stats)

    # leap year: February has 29 days
    @test stats.totrain2[2] == sum(obs.rain[(obs.year .== 2004) .& (obs.month .== 2)])

    # a missing day, a repeated day, an impossible date
    @test_throws r"year 2004 has 365 of its 366 days \(the first missing is 2004-02-29\)" M.weather_stats(
        obs[.!((obs.year .== 2004) .& (obs.month .== 2) .& (obs.day .== 29)), :])
    @test_throws r"2003-01-05 appears more than once" M.weather_stats(vcat(obs, obs[5:5, :]))
    bad = copy(obs)
    bad.day[40] = 30      # 9 February becomes 30 February
    @test_throws r"data row 40: 2003-2-30 is not a date" M.weather_stats(bad)
    @test_throws r"needs a column `day`" M.weather_stats(select(obs, Not(:day)))

    # missing or non-finite values
    bad = copy(obs)
    bad.tmax[10] = NaN
    @test_throws r"`tmax` has no valid value on 2003-01-10" M.weather_stats(bad)
    bad = allowmissing(obs)
    bad.rain[400] = missing
    @test_throws r"`rain` has no valid value on 2004-02-04" M.weather_stats(bad)
end


@testset "statistics of observed weather: constant month" begin
    obs = test_obs(2003)
    obs.tmax[obs.month .== 1] .= 30.0
    obs.wind[obs.month .== 2] .= 1.2
    stats = M.weather_stats(obs)
    @test stats.mean_tmax1[1] == 30.0
    @test stats.sd_tmax1[1] == 0.01        # sd is floored at 0.01
    @test stats.skew_tmax1[1] == 0.0       # not NaN
    @test stats.rlag_tmax1[1] == 0.0
    @test M.month_stats(fill(30.0, 31)).skew == 0.0

    # and weather can be generated from them
    nt = generate(stats)
    @test all(isfinite, nt.tmax[1].values)
    @test all(isfinite, nt.wind[1].values)
end


@testset "stats file checks" begin
    stats = M.weather_stats(test_obs(2003:2004))
    @test isnothing(M.check_stats(stats))

    function broken(col, v)
        s = copy(stats)
        s[!, col] = fill(v, nrow(s))
        s
    end
    @test_throws r"rlag_tmax3 of 2003 is 1.2; it must be between -1 and 1" M.check_stats(broken("rlag_tmax3", 1.2))
    @test_throws r"sd_tmin1 of 2003 is NaN; it must be above 0" M.check_stats(broken("sd_tmin1", NaN))
    @test_throws r"skew_tmin1 of 2003 is Inf" M.check_stats(broken("skew_tmin1", Inf))
    @test_throws r"mean_wind5 of 2003 is 0.0; it must be above 0" M.check_stats(broken("mean_wind5", 0.0))
    @test_throws r"totrain4 of 2003 is -1.0; it must be 0 or more" M.check_stats(broken("totrain4", -1.0))
    @test_throws r"pww2 of 2003 is 1.5; it must be between 0 and 1" M.check_stats(broken("pww2", 1.5))
    @test_throws r"column pwd7 is missing" M.check_stats(select(stats, Not(:pwd7)))
    @test_throws r"a year more than once" M.check_stats(vcat(stats, stats[1:1, :]))
    @test_throws r"has 4 problem" M.check_stats(broken("pww2", -0.5) |> s -> (s.pwd3 .= 2.0; s))

    # generation refuses an invalid stats file
    @test_throws r"rlag_tmax3" M.generate_mets(broken("rlag_tmax3", 1.2); verbose=false)
end


@testset "generation" begin
    stats = M.weather_stats(test_obs(2003:2005))

    # the same seed repeats a run exactly
    a = M.collate_mets(generate(stats; seed=7))
    b = M.collate_mets(generate(stats; seed=7))
    @test isequal(a, b)

    # whole consecutive years, every day once, tmin below tmax
    @test a.year == Dates.year.(Date(2003):Day(1):Date(2005, 12, 31))
    @test nrow(a) == 365 * 2 + 366     # 2004 is a leap year
    @test all(a.tmin .< a.tmax)
    @test all(a.wind .>= 0.1)
    @test all(a.rain .>= 0.0)

    # one fit score per month, for every variable and year
    nt = generate(stats)
    @test keys(nt) == (:tmin, :tmax, :wind, :rain)
    for mets ∈ values(nt), m ∈ mets
        @test length(m.scores) == 12
        @test all(s -> isfinite(s) && s >= 0, m.scores)
    end

    # a dry month stays dry, with no fit error
    dry = copy(stats)
    dry.totrain6 .= 0.0
    dry.pww6 .= 0.0
    dry.pwd6 .= 0.0
    nt = generate(dry)
    for m ∈ nt.rain
        june = M.month_ranges(m.year)[6]
        @test all(iszero, m.values[june])
        @test m.scores[6] == 0.0
    end

    # only some of the variables
    tmax_only = select(stats, :year, r"_tmax")
    nt = generate(tmax_only)
    @test keys(nt) == (:tmax,)
    @test names(M.collate_mets(nt)) == ["year", "month", "day", "doy", "tmax"]
end


@testset "fit diagnostics" begin
    # a month that cannot fit: the best score is returned after maxrun attempts
    data = zeros(31)
    nrun = Ref(0)
    s = M.autoregress_month!(data, 1:31, 0.0, 0.0, 0.0, e -> (nrun[] += 1; randn!(e)),
                             _ -> 5.0; maxrun=3)
    @test s == 5.0
    @test nrun[] == 3

    # a month that fits ends the search
    nrun[] = 0
    s = M.autoregress_month!(data, 1:31, 0.0, 0.0, 0.0, e -> (nrun[] += 1; randn!(e)),
                             _ -> 0.5)
    @test s == 0.5
    @test nrun[] == 1

    # months above 1 are listed, worst first
    met(year, scores) = M.Met{M.Temp}(year=year, scores=scores)
    nt = (tmax=[met(2003, [0.5; 1.4; fill(0.9, 10)]), met(2004, [3.0; fill(0.1, 11)])],)
    mf = M.misfits(nt)
    @test mf.year == [2004, 2003]
    @test mf.month == [1, 2]
    @test mf.score == [3.0, 1.4]
    txt = sprint(io -> M.print_misfits(io, nt))
    @test occursin("2 of 24 generated months missed the fit tolerance", txt)
    @test occursin("Tmax  2004 Jan   3.00", txt)
    ok = (tmax=[met(2003, fill(0.5, 12))],)
    @test occursin("All 12 generated months are within", sprint(io -> M.print_misfits(io, ok)))

    # Tmin >= Tmax is repaired and counted
    tn = [M.Met{M.Temp}(year=2003, values=[20.0, 25.0, 22.0])]
    tx = [M.Met{M.Temp}(year=2003, values=[30.0, 24.0, 22.0])]
    @test M.fix_tmin_tmax!(tn, tx) == 2
    @test tn[1].values == [20.0, 24.0, 22.0]
    @test tx[1].values == [30.0, 25.0, 22.1]
end


@testset "check_fit coverage" begin
    spec = M.weather_stats(test_obs(2003:2005))
    sim = M.weather_stats(M.collate_mets(generate(spec)))

    pairs, gaps = M.fit_pairs(spec, sim)
    @test isempty(gaps.years) && isempty(gaps.stats)
    @test sort(unique(pairs.year)) == [2003, 2004, 2005]

    # a missing year and a missing variable are reported, not hidden
    short = select(sim[sim.year .!= 2005, :], Not(r"_wind"))
    pairs, gaps = M.fit_pairs(spec, short)
    @test gaps.years == [2005]
    @test gaps.stats == ["mean_wind", "sd_wind", "rlag_wind"]
    @test 2005 ∉ pairs.year
    txt = sprint(io -> M.print_fit(io, "test", M.fit_summary(pairs), gaps))
    @test occursin("NOT COMPARED: the simulated weather lacks year(s) 2005", txt)
    @test occursin("NOT COMPARED: the simulated weather lacks wind, mean; wind, sd; wind, rlag", txt)

    @test_throws r"simulated weather has a year more than once" M.fit_pairs(spec, vcat(sim, sim[1:1, :]))
    @test_throws r"no year or statistic" M.fit_pairs(spec, sim[sim.year .== 1999, :])

    # tolerances, in each statistic's own units
    @test M.stat_tol("mean_tmax", 32.0) == M.TEMP_TOL.mean
    @test M.stat_tol("sd_tmin", 0.8) ≈ M.TEMP_TOL.sd * 0.8
    @test M.stat_tol("totrain", 200.0) ≈ 5.0
    @test M.stat_tol("pww", 0.5) ≈ 0.025
    @test M.stat_tol("pwd", 0.0) ≈ 0.0005

    summ = M.fit_summary(M.fit_pairs(spec, sim)[1])
    @test all(0 .<= summ.n_out .<= summ.n)
end


@testset "whole run: Serdang" begin
    src = joinpath(@__DIR__, "..", "data", "Serdang")
    mktempdir() do folder
        mkpath(joinpath(folder, "Serdang"))
        for f ∈ ("Serdang-obs.csv", "Serdang-stats.csv")
            cp(joinpath(src, f), joinpath(folder, "Serdang", f))
        end
        df = redirect_stdout(devnull) do
            generate_weather("Serdang"; folder=folder, seed=1)
        end
        @test nrow(df) == 365 * 4 + 366
        summ = redirect_stdout(devnull) do
            check_fit("Serdang"; folder=folder, show=false)
        end
        @test :n_out ∈ propertynames(summ)
        @test isfile(joinpath(folder, "Serdang", "Serdang-fit.txt"))
    end
end

end
