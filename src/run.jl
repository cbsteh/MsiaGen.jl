
"""
    generate_weather(site; folder="data", seed=-1, from_stats=false, verbose=false,
                     wetdry_tmax=-0.6, io=stdout)

Generate daily weather for `site`, whose files are in `folder/site/`:

- `<site>-obs.csv`: observed daily weather (input)
- `<site>-stats.csv`: monthly and annual statistics of the observed weather
- `<site>-sim.csv`: simulated daily weather (output)

The statistics file is rebuilt from the observed file, unless `from_stats=true`,
in which case an existing `<site>-stats.csv` is used (no observed file needed).
A negative `seed` picks a random seed. Rain is generated first, and Tmax is
`wetdry_tmax` (°C) warmer on wet days than on dry days of the same month
(negative: cooler; 0 for no difference, as in the published version), with
each month's Tmax statistics unchanged. Prints the months (if any) whose
generated weather missed the fit tolerance, to `io`. Returns the simulated weather.
"""
function generate_weather(
    site::AbstractString;
    folder::AbstractString="data",
    seed::Int=-1,
    from_stats::Bool=false,
    verbose::Bool=false,
    wetdry_tmax::Real=WETDRY_TMAX,
    io::IO=stdout
)
    path = joinpath(folder, site)
    obs_file = joinpath(path, "$(site)-obs.csv")
    stats_file = joinpath(path, "$(site)-stats.csv")
    sim_file = joinpath(path, "$(site)-sim.csv")

    from_stats || create_data_file(obs_file, basename(stats_file))

    seednum = (seed < 0) ? rand(1:typemax(Int)) : seed
    Random.seed!(seednum)
    println(io, ">>> Generating for $(site). Using seed no. $(seednum).")

    stats = csv2df(stats_file)
    nt = generate_mets(stats.df; verbose=verbose, wetdry_tmax=wetdry_tmax, io=io)
    print_misfits(io, nt)
    df = collate_mets(nt)
    CSV.write(sim_file, df)

    println(io, "...written to $(sim_file)")
    df
end
