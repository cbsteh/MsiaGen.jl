
"""
    generate_weather(site; folder="data", seed=-1, from_stats=false, verbose=false)

Generate daily weather for `site`, whose files are in `folder/site/`:

- `<site>-obs.csv`: observed daily weather (input)
- `<site>-stats.csv`: monthly and annual statistics of the observed weather
- `<site>-sim.csv`: simulated daily weather (output)

The statistics file is rebuilt from the observed file, unless `from_stats=true`,
in which case an existing `<site>-stats.csv` is used (no observed file needed).
A negative `seed` picks a random seed. Returns the simulated weather.
"""
function generate_weather(
    site::AbstractString;
    folder::AbstractString="data",
    seed::Int=-1,
    from_stats::Bool=false,
    verbose::Bool=false
)
    path = joinpath(folder, site)
    obs_file = joinpath(path, "$(site)-obs.csv")
    stats_file = joinpath(path, "$(site)-stats.csv")
    sim_file = joinpath(path, "$(site)-sim.csv")

    from_stats || create_data_file(obs_file, basename(stats_file))

    seednum = (seed < 0) ? rand(1:typemax(Int)) : seed
    Random.seed!(seednum)
    println(">>> Generating for $(site). Using seed no. $(seednum).")

    stats = csv2df(stats_file)
    df = collate_mets(generate_mets(stats.df; verbose=verbose))
    CSV.write(sim_file, df)

    println("...written to $(sim_file)")
    df
end
