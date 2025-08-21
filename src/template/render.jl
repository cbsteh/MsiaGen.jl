
function replace_words(str::AbstractString, dict::AbstractDict)
    for (key, value) in dict
        key_with_braces = r"{{\s*" * key * r"\s*}}"
        str = replace(str, key_with_braces => value)
    end
    str
end


function encode64(fig::Figure)
    io = IOBuffer()
    iob64_encode = Base64EncodePipe(io)
    show(iob64_encode, MIME("image/png"), fig)
    close(iob64_encode)
    String(take!(io))
end


function fig2html(fig::Figure)
    data64 = encode64(fig)
    img = "<img src=\"data:image/png;base64,$(data64)\" alt=\"\"/>"
    cap = "<figcaption class=\"hr-lines\">\u2609 \u2609</figcaption>"
    "<figure>$(img)$(cap)</figure>"
end


function fig2html(fig_lst::AbstractVector)
    html = ""
    for fig ∈ fig_lst
        if !isnothing(fig)
            html *= fig2html(fig)
        end
    end

    if isempty(html)
        html = "<p>No plots</p>"
    end

    html
end


function render2html(
    tmp_file::AbstractString,
    html_file::AbstractString,
    dict::AbstractDict
)
    open(tmp_file, "r") do fin
        content = read(fin, String)
    end

    content = replace_words(content, dict)

    open(html_file, "w") do fout
        write(fout, content)
    end
end


function report(
    obs_file::AbstractString,
    data_file::AbstractString,
    gen_file::AbstractString;
    verbose::Bool=true,
    create_gen_file::Bool=true
)
    data = csv2df(data_file)
    nt = generate_mets(data.df; verbose=verbose)
    est = collate_mets(nt)

    create_gen_file && CSV.write(gen_file, est)

    add_solar_radiation!(est, deg2rad(data.lat))

    obs = csv2df(obs_file)
    stat = collate_stats(nt; transpose=true)
    pfit = gof(nt; ndecimals=2)
    obs.df, est, stat, pfit
end


function fit_dist_metrics(obs, est, site, k1, k2, k3, k4)
    function ADtest(o, e, param)
        gobs = groupby(o, :year)
        gest = groupby(e, :year)
        ps = Float64[]
        for (i, j) in zip(eachindex(gest), eachindex(gobs))
            x = getproperty(gest[i], param)
            y = getproperty(gobs[j], param)
            p = pvalue(KSampleADTest(x, y))
            push!(ps, p)
        end
        ps
    end

    function KStest(o, e, param)
        gobs = groupby(o, :year)
        gest = groupby(e, :year)
        ps = Float64[]
        for (i, j) in zip(eachindex(gest), eachindex(gobs))
            x = getproperty(gest[i], param)
            y = getproperty(gobs[j], param)
            p = pvalue(ApproximateTwoSampleKSTest(x, y))
            push!(ps, p)
        end
        ps
    end

    results = DataFrame()
    params = [
        (:tmin, k1, false),
        (:tmax, k2, false),
        (:wind, k3, false),
        (:rain, k4, true)
    ]

    for (col, flag, rainfilter) in params
        isnothing(flag) && continue

        o = obs
        e = est
        if rainfilter
            o = obs[obs.rain .> 0, :]
            e = est[est.rain .> 0, :]
        end

        ad_res = DataFrame(AD = ADtest(o, e, col), KS = KStest(o, e, col))
        insertcols!(ad_res, 1, :PARAM=>string(col))
        results = vcat(results, ad_res)
    end

    insertcols!(results, 1, :SITE=>site)
    results
end


function fit_monthly_metrics(obs, est, site)
    gobs = groupby(obs, :year)
    gest = groupby(est, :year)

    metrics = [GoF.NMAE, GoF.NMBE, GoF.KGE, GoF.dr]
    str_metrics = Symbol.(metrics)

    has = exists_met(est)
    params = [
        (has.tmin, mean, :tmin),
        (has.tmax, mean, :tmax),
        (has.wind, mean, :wind),
        (has.rain, sum,  :rain)
    ]

    sbm = summarize_by_month
    results = DataFrame()

    for (enabled, agg, param) in params
        !enabled && continue
        for (i, j) in zip(eachindex(gest), eachindex(gobs))
            x = getproperty(gest[i], param)
            y = getproperty(gobs[j], param)
            year = gobs[j].year[1]

            xmth = sbm(year, Vector(x), agg)
            ymth = sbm(year, Vector(y), agg)
            vals = (g -> g(ymth, xmth)).(metrics)

            nt = (; (k=>[v] for (k, v) in zip(str_metrics, vals))...)
            ad = DataFrame(nt)
            insertcols!.(Ref(ad), [1,2], [:SITE=>site, :PARAM=>string(param)])
            results = vcat(results, ad)
        end
    end

    results
end



function report_to_screen(
    obs_file::AbstractString,
    data_file::AbstractString,
    gen_file::AbstractString,
    site::AbstractString;
    verbose::Bool=true,
    create_gen_file::Bool=true
)
    obs, est, stat, pfit = report(obs_file, data_file, gen_file;
                                  verbose=verbose, create_gen_file=create_gen_file)
    plot_monthly(obs, est; plot=true, site=site)
    k1, k2, k3, k4 = plot_kde(obs, est; plot=true, site=site)
    plot_r(obs, est; plot=true, site=site)
    plot_daily(obs, est; plot=true, site=site)
    plot_solar(est; plot=true, site=site)

    mfit = fit_monthly_metrics(obs, est, site)
    ad = fit_dist_metrics(obs, est, site, k1, k2, k3, k4)

    (; obs=obs, est=est, stat=stat, pfit=pfit, mfit=mfit, ad=ad)
end


function report_to_html(
    obs_file::AbstractString,
    data_file::AbstractString,
    gen_file::AbstractString,
    site::AbstractString;
    verbose::Bool=true,
    create_gen_file::Bool=true
)
    obs, est, stat, pfit = report(obs_file, data_file, gen_file;
                                  verbose=verbose, create_gen_file=create_gen_file)

    plot = false
    m1, m2, m3, m4 = plot_monthly(obs, est; site=site, plot=plot)
    k1, k2, k3, k4 = plot_kde(obs, est; site=site, plot=plot)
    r1 = plot_r(obs, est; site=site, plot=plot)
    d1, d2, d3 = plot_daily(obs, est; site=site, plot=plot)

    figs_monthly = fig2html([m1, m2, m3, m4])
    figs_kde = fig2html([k1, k2, k3, k4])
    figs_r = fig2html([r1])
    figs_daily = fig2html([d1, d2, d3])

    fit_table = pretty_table(String, pfit; backend=Val(:html),
                             show_subheader=false, table_class="minimalistBlack")
    stat_table = pretty_table(String, stat; backend=Val(:html),
                              show_subheader=false, table_class="minimalistBlack")
    loc = replace(site, "-" => " ", "_" => " ")

    dict = Dict(
        "title" => loc,
        "figs_monthly" => figs_monthly,
        "figs_kde" => figs_kde,
        "figs_r" => figs_r,
        "figs_daily" => figs_daily,
        "fit_table" => fit_table,
        "stat_table" => stat_table
    )

    html_file = dirname(data_file) * "/$(site).html"
    render2html("src/template/report.html", html_file, dict)

    annfit = fit_monthly_metrics(obs, est, site)
    distfit = fit_dist_metrics(obs, est, site, k1, k2, k3, k4)

    (; obs=obs, est=est, stat=stat, pfit=pfit, annfit=annfit, distfit=distfit)
end


function batch_report_to_html(
    basedir::AbstractString;
    has_data_file::Bool=true,
    create_gen_file::Bool=true,
    res_annual_fit::AbstractString="results/fits/_annfit.csv",
    res_dist_fit::AbstractString="results/fits/_distfit.csv"
)
    ops = Vector{String}[]

    for (root, dirs, files) in walkdir(basedir)
        isempty(files) && continue

        site = splitpath(root)[end]
        obs_file = joinpath(root, "$(site).csv")
        gen_file = joinpath(root, "wthr.csv")
        if isfile(obs_file)
            data_file = has_data_file ?
                joinpath(root, "data.csv") :
                create_data_file(obs_file, "data.csv")
            isfile(data_file) && push!(ops, [obs_file, data_file, gen_file, site])
        end
    end

    length(ops) == 0 && return

    println("From \"$(basedir)\", preparing a report for:")
    rm(res_annual_fit; force=true)
    rm(res_dist_fit; force=true)

    for (cnt, op) ∈ enumerate(ops)
        println("\t$(op[4])")
        bheader = (cnt==1)
        out = report_to_html(op[1], op[2], op[3], op[4];
                             verbose=false, create_gen_file=create_gen_file)
        CSV.write(res_annual_fit, out.annfit, append=true, header=bheader)
        CSV.write(res_dist_fit, out.distfit, append=true, header=bheader)
    end

    println("...done.")
end
