module MsiaGen

using CSV
using DataFrames
using Dates
using Distributions
using Parameters
using Printf
using Random
using SpecialFunctions
using Statistics
using StatsBase

include("met.jl")
include("utils.jl")
include("collate.jl")
include("gentemp.jl")
include("genwind.jl")
include("genrain.jl")
include("autoreg.jl")
include("input.jl")
include("run.jl")
include("plotwthr.jl")

using .Plotting

include("checkfit.jl")

# utils.jl
export csv2df
# collate.jl
export collate_mets, generate_mets
# gentemp.jl
export create_temp
# genwind.jl
export create_wind
# genrain.jl
export create_rain
# gentemp.jl, genwind.jl, genrain.jl
export generate!
# input.jl
export create_data_file
# run.jl
export generate_weather
# plotwthr.jl
export plot_weather
# checkfit.jl
export check_fit


# Run the whole pipeline once on a small made-up site while the package
# precompiles, so the first real run does not wait for compilation
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    days = 1:365
    obs = DataFrame(year=2001, month=Dates.month.(Date(2001) .+ Day.(days .- 1)),
                    day=Dates.day.(Date(2001) .+ Day.(days .- 1)),
                    tmin=23.0 .+ sin.(days ./ 7) .+ 0.3 .* cos.(days .* 1.3),
                    tmax=32.0 .+ sin.(days ./ 5) .+ 0.5 .* cos.(days .* 1.7),
                    wind=1.5 .+ 0.3 .* sin.(days ./ 3) .+ 0.2 .* cos.(days .* 2.1),
                    rain=[i % 3 == 0 ? 2.0 + i % 11 : 0.0 for i ∈ days])
    @compile_workload begin
        mktempdir() do folder
            site = "site"
            mkpath(joinpath(folder, site))
            open(joinpath(folder, site, "$(site)-obs.csv"), "w") do io
                println(io, 3.0)
                CSV.write(io, obs; append=true, writeheader=true)
            end
            redirect_stdout(devnull) do
                generate_weather(site; folder=folder, seed=1)
                plot_weather(site; folder=folder, show=false)
                check_fit(site; folder=folder, show=false)
            end
        end
    end
end

end # module MsiaGen
