
"""
    check_links(site; folder="data")

Compare how the weather variables move together in the observed
(`<site>-obs.csv`) and simulated (`<site>-sim.csv`) weather. MsiaGen
generates each variable on its own, so the simulated weather has these
links only by chance. For each month, and all months together:

- Tmin-Tmax r: correlation of the daily tmin and tmax
- sd(Tmax-Tmin): sd of the daily temperature range (°C)
- Tmax wet-dry, Tmin wet-dry: mean on wet days minus mean on dry days (°C);
  a wet day has rain above 0

Each measure is worked out within each month of each year (so differences
between years do not count), then averaged over the years. Prints the
observed and simulated values side by side and writes them to the plain
text file `<site>-links.txt`. Returns the table.
"""
function check_links(site::AbstractString; folder::AbstractString="data")
    path = joinpath(folder, site)
    obs = weather_links(csv2df(joinpath(path, "$(site)-obs.csv")).df)
    sim = weather_links(DataFrame(CSV.File(joinpath(path, "$(site)-sim.csv"))))
    tbl = DataFrame(month=obs.month)
    for k ∈ keys(LINKS)
        tbl[!, "obs_$(k)"] = obs[!, k]
        tbl[!, "sim_$(k)"] = sim[!, k]
    end
    report = sprint(io -> print_links(io, site, tbl))
    print(report)
    txt = joinpath(path, "$(site)-links.txt")
    write(txt, report)
    println("\nTable saved as $(txt)")
    tbl
end


# The measures of check_links: column name and heading
const LINKS = (r_tmin_tmax="Tmin-Tmax r", sd_dtr="sd(Tmax-Tmin)",
               wet_tmax="Tmax wet-dry", wet_tmin="Tmin wet-dry")


# The measures of one month of one year of daily weather `g`
function month_links(g)
    wet = g.rain .> 0
    gap(x) = any(wet) && !all(wet) ? mean(x[wet]) - mean(x[.!wet]) : NaN
    (r_tmin_tmax=cor(g.tmin, g.tmax), sd_dtr=std(g.tmax .- g.tmin),
     wet_tmax=gap(g.tmax), wet_tmin=gap(g.tmin))
end


# The measures of the daily weather `df` for each month (1-12) and all
# months (0), averaged over the years; a month-year whose measure cannot be
# worked out (no wet or no dry day, constant temperature) is left out
function weather_links(df::AbstractDataFrame)
    weather_dates(df)
    for c ∈ ("tmin", "tmax", "rain")
        c ∈ names(df) || error("check_links needs the column `$(c)`")
    end
    per = combine(groupby(df, [:year, :month]),
                  AsTable([:tmin, :tmax, :rain]) => month_links => AsTable)
    avg(x) = (y = filter(!isnan, x); isempty(y) ? NaN : mean(y))
    row(m, d) = (month=m, (k => avg(d[!, k]) for k ∈ keys(LINKS))...)
    DataFrame([[row(m, per[per.month .== m, :]) for m ∈ 1:12]; row(0, per)])
end


function print_links(io::IO, site, tbl)
    println(io, "\nLinks between the weather variables, $(site): observed (obs) " *
                "and simulated (sim)")
    println(io, "Each measure is worked out within each month of each year, then " *
                "averaged over the years.\n")
    @printf(io, "%-5s", "Month")
    foreach(h -> @printf(io, " %15s", h), values(LINKS))
    @printf(io, "\n%-5s", "")
    foreach(_ -> @printf(io, " %7s %7s", "obs", "sim"), LINKS)
    println(io)
    for r ∈ eachrow(tbl)
        @printf(io, "%-5s", r.month == 0 ? "All" : MONTHS[r.month])
        for k ∈ keys(LINKS)
            @printf(io, " %+7.2f %+7.2f", r["obs_$(k)"], r["sim_$(k)"])
        end
        println(io)
    end
    println(io, "\nsd(Tmax-Tmin), Tmax wet-dry and Tmin wet-dry are in °C; a wet " *
                "day has rain above 0.")
end
