
"""
    check_fit(site; folder="data", show=true)

Check how well the simulated weather (`<site>-sim.csv`) matches the specified
statistics (`<site>-stats.csv`). The simulated weather is measured with the
same method that builds the stats file (`create_data_file`), and every
monthly statistic is compared with its specified value, for every year.

Prints, for each statistic, the mean absolute error, the bias (mean of
simulated minus specified), the worst error (with its month and year), all
in the statistic's own units, and how many values are outside the
generator's fit tolerance. Years or statistics of the stats file that the
simulated weather lacks are listed above the table. Writes the same report to the plain
text file `<site>-fit.txt`. Saves a chart of specified against
simulated values, one panel per statistic, as `<site>-fit.png`, displayed
if `show=true`. Returns the summary table.
"""
function check_fit(site::AbstractString; folder::AbstractString="data", show::Bool=true)
    path = joinpath(folder, site)
    spec = csv2df(joinpath(path, "$(site)-stats.csv"))
    sim = weather_stats(DataFrame(CSV.File(joinpath(path, "$(site)-sim.csv"))))

    pairs, gaps = fit_pairs(spec.df, sim)
    summ = fit_summary(pairs)
    report = sprint(io -> print_fit(io, site, summ, gaps))
    print(report)
    txt = joinpath(path, "$(site)-fit.txt")
    write(txt, report)

    pairs.label = stat_label.(pairs.stat, pairs.month)
    pairs.unit = stat_unit.(pairs.stat)
    pairs.min_span = min_span.(pairs.stat, pairs.spec)
    fname = joinpath(path, "$(site)-fit.png")
    Plotting.plot_fit(site, pairs, fname; show=show)
    println("\nTable saved as $(txt); chart saved as $(fname)")
    summ
end


# One row per statistic, month and year: specified and simulated values.
# Monthly values only (months 1-12), plus the annual rain total (month 0).
# Also returns what the simulated weather lacks: years and statistics of
# the stats file it does not have.
function fit_pairs(spec::AbstractDataFrame, sim::AbstractDataFrame)
    allunique(spec.year) || error("the stats file has a year more than once")
    allunique(sim.year) || error("the simulated weather has a year more than once")
    row_of(df) = Dict(y => i for (i, y) ∈ enumerate(df.year))
    spec_row, sim_row = row_of(spec), row_of(sim)
    years = sort(intersect(spec.year, sim.year))
    rows = NamedTuple[]
    lacking = String[]
    for col ∈ names(spec)
        m = match(r"^(.*?)(\d+)$", col)
        isnothing(m) && continue
        stat, month = m[1], parse(Int, m[2])
        (month == 0 && stat != "totrain") && continue
        if col ∉ names(sim)
            stat ∈ lacking || push!(lacking, stat)
            continue
        end
        a, b = spec[!, col], sim[!, col]
        for y ∈ years
            push!(rows, (stat=stat, month=month, year=y,
                         spec=Float64(a[spec_row[y]]), sim=Float64(b[sim_row[y]])))
        end
    end
    isempty(rows) && error("the simulated weather has no year or statistic of the stats file")
    gaps = (years=sort(setdiff(spec.year, sim.year)), stats=lacking)
    DataFrame(rows), gaps
end


# Readable name and unit of a statistic's column prefix
function stat_label(stat::AbstractString, month::Int=1)
    stat == "totrain" && return month == 0 ? "rain, annual total" : "rain, monthly total"
    stat ∈ ("pww", "pwd") && return "rain, $(stat)"
    p, v = split(stat, "_")
    "$(v), $(p)"
end


# Smallest axis range of a statistic's panel in the fit chart: 10 times the
# generator's fit tolerance, so a statistic whose specified value never
# changes shows as a small cluster on the 1:1 line, not a tall column
function min_span(stat::AbstractString, value::Real)
    k = 10
    rel(tol) = k * tol * abs(value)
    startswith(stat, "mean_t") ? k * TEMP_TOL.mean :
    stat == "mean_wind" ? k * WIND_TOL.mean :
    startswith(stat, "sd_t") ? rel(TEMP_TOL.sd) :
    stat == "sd_wind" ? rel(WIND_TOL.sd) :
    startswith(stat, "rlag_t") ? k * TEMP_TOL.rlag :
    stat == "rlag_wind" ? k * WIND_TOL.rlag :
    startswith(stat, "skew") ? k * TEMP_TOL.skew :
    stat == "totrain" ? rel(0.025) :          # rain fit tolerance: 2.5% of the total
    k * 0.05                                  # pww, pwd
end


# The generator's fit tolerance of a statistic whose specified value is
# `value`, in the statistic's own units (gentemp.jl, genwind.jl, genrain.jl)
function stat_tol(stat::AbstractString, value::Real)
    rel(pct) = pct / 100 * max(abs(value), 0.01)
    startswith(stat, "mean_t") ? TEMP_TOL.mean :
    stat == "mean_wind" ? WIND_TOL.mean :
    startswith(stat, "sd_t") ? TEMP_TOL.sd * value :
    stat == "sd_wind" ? WIND_TOL.sd * value :
    startswith(stat, "rlag_t") ? TEMP_TOL.rlag :
    stat == "rlag_wind" ? WIND_TOL.rlag :
    startswith(stat, "skew") ? TEMP_TOL.skew :
    stat == "totrain" ? rel(RAIN_TOL.totrain) :
    rel(RAIN_TOL.pw)                          # pww, pwd
end


function stat_unit(stat::AbstractString)
    stat == "totrain" && return "mm"
    startswith(stat, "mean_t") || startswith(stat, "sd_t") ? "°C" :
        stat ∈ ("mean_wind", "sd_wind") ? "m/s" : ""
end


# Mean absolute error, bias, worst error and the number of values outside
# the fit tolerance, of each statistic
function fit_summary(pairs::AbstractDataFrame)
    pairs = transform(pairs, [:sim, :spec] => ((s, p) -> s .- p) => :err,
                      [:stat, :spec] => ByRow(stat_tol) => :tol,
                      [:stat, :month] => ByRow((s, m) -> s == "totrain" && m == 0 ?
                                                         "totrain0" : s) => :group)
    combine(groupby(pairs, :group; sort=false)) do g
        i = argmax(abs.(g.err))
        pct = g.stat[1] == "totrain" ?
              100 * mean(abs.(g.err) ./ max.(g.spec, 1.0)) : missing
        (label=stat_label(g.stat[1], g.month[1]), unit=stat_unit(g.stat[1]),
         n=nrow(g), mae=mean(abs.(g.err)), bias=mean(g.err), mae_pct=pct,
         n_out=count(abs.(g.err) .> g.tol),
         worst=g.err[i], worst_when=(g.month[i] == 0 ? "" : MONTHS[g.month[i]] * " ") *
                                    string(g.year[i]))
    end
end


function print_fit(io::IO, site, summ, gaps=(years=Int[], stats=String[]))
    println(io, "\nFit of simulated to specified statistics, $(site) " *
                "(error = simulated - specified)\n")
    isempty(gaps.years) ||
        println(io, "NOT COMPARED: the simulated weather lacks year(s) " *
                    join(gaps.years, ", ") * " of the stats file\n")
    isempty(gaps.stats) ||
        println(io, "NOT COMPARED: the simulated weather lacks " *
                    join(stat_label.(gaps.stats), "; ") * "\n")
    @printf(io, "%-22s %15s %10s %22s %12s\n", "Statistic", "Mean |error|", "Bias",
            "Worst error", "Outside tol")
    for r ∈ eachrow(summ)
        u = isempty(r.unit) ? "" : " " * r.unit
        mae = @sprintf("%.3g%s", r.mae, u) *
              (ismissing(r.mae_pct) ? "" : @sprintf(" (%.1f%%)", r.mae_pct))
        @printf(io, "%-22s %15s %+10.3g %22s %12s\n", r.label, mae, r.bias,
                @sprintf("%+.3g%s (%s)", r.worst, u, r.worst_when), "$(r.n_out)/$(r.n)")
    end
    println(io, "\nOutside tol: values further from the specified value than the " *
                "generator's fit tolerance.")
end
