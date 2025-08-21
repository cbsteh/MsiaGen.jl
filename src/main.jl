using MsiaGen
using Random


function generate_weather(
    folder::AbstractString,
    site::AbstractString;
    seed::Int=-1,
    has_data_file::Bool=:false,
    verbose::Bool=false
)
    path = joinpath(folder, site)
    obs_file = joinpath(path, "$(site).csv")
    data_file = has_data_file ?
                    joinpath(path, "data.csv") :
                    create_data_file(obs_file, "data.csv")
    gen_file = joinpath(path, "wthr.csv")

    seednum = (seed < 0) ? rand(1:typemax(Int)) : seed
    Random.seed!(seednum)

    println(">>> Generating for $(site). Using seed no. $seednum.")

    obs, est, stat, pfit = report(obs_file, data_file, gen_file; verbose=verbose)

    plot_monthly(obs, est; plot=true, site=site)
    k1, k2, k3, k4 = plot_kde(obs, est; plot=true, site=site)
    plot_r(obs, est; plot=true, site=site)

    mfit = fit_monthly_metrics(obs, est, site)
    ad = fit_dist_metrics(obs, est, site, k1, k2, k3, k4)

    println("...done.")

    (; obs=obs, est=est, stat=stat, pfit=pfit, mfit=mfit, ad=ad)
end



folder = "data"
site = "Serdang"
has_data_file = false
seed = -1
verbose = false

res = generate_weather(
    folder, site;
    seed=seed,
    has_data_file=has_data_file,
    verbose=verbose
)

;
