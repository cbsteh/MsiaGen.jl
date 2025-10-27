include("nicescale.jl")


function num_rows_cols(N, maxcols)
    n, r = divrem(N, maxcols)
    ncols = (n > 0) ? maxcols : r
    nrows = n + (r>0)
    nrows, ncols
end


function nextpos(row, col, maxcols)
    col += 1
    col = max(1, col % (maxcols+1))
    if col == 1
        row += 1
    end
    row, col
end


function add_legend_metrics!(ax, x, y, metrics, i)
    vals = round.((g -> g(y, x)).(metrics); digits=2)
    key = (i==1) ? join(string.(metrics), ", ") * "\n" : ""
    labels = key * join(vals, ", ")
    dummy = scatter!(ax, [x[1]], [y[1]], color=:transparent, markersize=1)
    axislegend(
        ax, [dummy], [labels],
        padding=(-10, 10, 1, 1),
        colgap=3,
        orientation=:horizontal,
        position=:rb,
        backgroundcolor=:transparent,
    )
end


function make_equal_scale!(ax, x, y)
    xy = [x; y]
    mn = minimum(xy)
    mx = maximum(xy)
    Δ = (mx - mn) * 0.1     # leave a 10% margin
    mn -= Δ
    mx += Δ
    limits!(ax, (mn, mx), (mn, mx))
end


function create_fig(N, maxcols, title)
    nrows, ncols = num_rows_cols(N, maxcols)
    dpi = 300
    fontsize = 18
    chtsize = 2    # inches
    sz_px = (chtsize * ncols * dpi, chtsize * nrows * dpi)
    fig = Figure(; size=sz_px, font="Arial", fontsize=fontsize)
    Label(fig[0, 1:ncols], title, font=:bold)   # chart title

    fig, nrows, ncols
end


function plot_daily(obs, est, param::Symbol, title::AbstractString, plot::Bool=true)
    gobs = groupby(obs, :year)
    gest = groupby(est, :year)

    N = length(gobs)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)
    metrics = [GoF.NMAE, GoF.NMBE, GoF.KGE, GoF.dr]

    axs = []
    row = col = 1
    for i ∈ 1:N
        x = getproperty(gest[i], param)
        y = getproperty(gobs[i], param)
        year = gobs[i].year[1]

        push!(axs, Axis(
            fig[row, col],
            xgridvisible=false, ygridvisible=false,
            xlabel="gen.", ylabel="obs.",
            title=string(year))
        )

        scatter!(axs[end], x, y, color=:black, markersize=7)
        ablines!(axs[end], 0, 1, color=:red, linestyle=:dash, linewidth=3)
        add_legend_metrics!(axs[end], x, y, metrics, i)
        make_equal_scale!(axs[end], x, y)

        row, col = nextpos(row, col, maxcols)
    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)
    fig
end


function plot_monthly(obs, est, fn::Function, param::Symbol, title::AbstractString,
                      plot::Bool=true)
    gobs = groupby(obs, :year)
    gest = groupby(est, :year)

    N = length(gobs)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)
    metrics = [GoF.NMAE, GoF.NMBE, GoF.KGE, GoF.dr]
    sbm = summarize_by_month

    axs = []
    row = col = 1
    for i ∈ 1:N
        x = getproperty(gest[i], param)
        y = getproperty(gobs[i], param)
        year = gobs[i].year[1]

        push!(axs, Axis(
            fig[row, col],
            xgridvisible=false, ygridvisible=false,
            xlabel="gen.", ylabel="obs.",
            title=string(year))
        )

        xmth = sbm(year, Vector(x), fn)
        ymth = sbm(year, Vector(y), fn)
        scatter!(axs[end], xmth, ymth, color=:black, markersize=12)

        ablines!(axs[end], 0, 1, color=:red, linestyle=:dash, linewidth=3)
        add_legend_metrics!(axs[end], xmth, ymth, metrics, i)
        make_equal_scale!(axs[end], xmth, ymth)

        row, col = nextpos(row, col, maxcols)
    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)
    fig
end


function psig(pval)
    sig = (pval < 0.01) ? "**" :
          (pval < 0.05) ? "*" : ""
    "p<$(round(pval; digits=3))$(sig)"
end


function plot_kde(obs, est, param::Symbol, title::AbstractString, plot::Bool=true)
    gobs = groupby(obs, :year)
    gest = groupby(est, :year)

    N = length(gobs)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)

    axs = []
    row = col = 1
    for i ∈ 1:N
        x = getproperty(gest[i], param)
        y = getproperty(gobs[i], param)
        year = gobs[i].year[1]

        p = pvalue(KSampleADTest(x, y))
        title = "$year ($(psig(p)))"

        push!(axs, Axis(
            fig[row, col],
            xgridvisible=false, ygridvisible=false,
            ylabel=string(param), title=title)
        )

        x1 = repeat([1], length(x))
        lmt = extrema(x)
        violin!(axs[end], x1, x; color=:teal, strokewidth=1, width=0.5,
                orientation=:vertical, datalimits=lmt, side=:left, show_median=true)

        x2 = repeat([1], length(y))
        lmt = extrema(y)
        violin!(axs[end], x2, y; color=:orange, strokewidth=1, width=0.5,
                orientation=:vertical, datalimits=lmt, side=:right, show_median=true)

        b1 = boxplot!(axs[end], x1 .- 0.5, x;
                      orientation=:vertical, width=0.1,
                      color = (:teal, 1), strokewidth=1,
                      markersize=7, outliercolor=:black)

        b2 = boxplot!(axs[end], x2 .+ 0.5, y;
                      orientation=:vertical, width=0.1,
                      color = (:orange, 1), strokewidth=1,
                      markersize=7, outliercolor=:black)

        axs[end].yticks = niceticks(vcat(x, y), 6)
        hidexdecorations!(axs[end])

        i == N && Legend(fig[row+1, :], [b1, b2], ["gen.", "obs."],
                         orientation=:horizontal, framevisible=false,
                         colgap=10, backgroundcolor=:transparent)
        row, col = nextpos(row, col, maxcols)

    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)

    fig
end


function linreg(x1, y1, x2, y2)
    df1 = DataFrame(x1=x1, y1=y1)
    df2 = DataFrame(x2=x2, y2=y2)
    model1 = lm(@formula(y1 ~ x1), df1)
    model2 = lm(@formula(y2 ~ x2), df2)
    c1, b1 = coef(model1)
    c2, b2 = coef(model2)
    X1 = [extrema(x1)...]
    X2 = [extrema(x2)...]
    Y1 = c1 .+ b1 .* X1
    Y2 = c2 .+ b2 .* X2

    sb1 = stderror(model1)[2]
    sb2 = stderror(model2)[2]
    t = (b1 - b2) / sqrt(sb1^2 + sb2^2)
    degf = size(df1, 1) + size(df2, 1) - 4
    p = 2 * ccdf(TDist(degf), abs(t))

    (; X1=X1, Y1=Y1, X2=X2, Y2=Y2, m1=b1, m2=b2, c1=c1, c2=c2, p=p, t=t, df=degf)
end


function plot_r(obs, est,
                params::AbstractVector{Symbol}, fn::AbstractVector{Function},
                title::AbstractString, plot::Bool=true)
    gobs = groupby(obs, :year)
    gest = groupby(est, :year)

    N = length(gobs)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)
    sbm = summarize_by_month

    axs = []
    row = col = 1
    for i ∈ 1:N
        x1 = getproperty(gest[i], params[1])
        y1 = getproperty(gest[i], params[2])
        x2 = getproperty(gobs[i], params[1])
        y2 = getproperty(gobs[i], params[2])
        year = gobs[i].year[1]

        x1mth = sbm(year, Vector(x1), fn[1])
        y1mth = sbm(year, Vector(y1), fn[2])
        x2mth = sbm(year, Vector(x2), fn[1])
        y2mth = sbm(year, Vector(y2), fn[2])

        reg = linreg(x1mth, y1mth, x2mth, y2mth)
        r1 = round(cor(x1mth, y1mth); digits=2)
        r2 = round(cor(x2mth, y2mth); digits=2)

        title = string(year)
        label1 = "gen. ($r1)"
        label2 = "obs. ($r2)"

        push!(axs, Axis(
            fig[row, col],
            xgridvisible=false, ygridvisible=false,
            xlabel=string(params[1]), ylabel=string(params[2]),
            title=title)
        )

        s1 = scatter!(axs[end], x1mth, y1mth, color=:red, markersize=12)
        lines!(axs[end], reg.X1, reg.Y1, linewidth=2, color=:red)

        s2 = scatter!(axs[end], x2mth, y2mth, color=:black, markersize=12, marker=:xcross)
        lines!(axs[end], reg.X2, reg.Y2, linewidth=2, color=:black)

        axislegend(axs[end], [s1, s2], [label1, label2], padding=(0, 2, 0, 0),
                   rowgap=-3, orientation=:vertical, position=:rt,
                   backgroundcolor=:transparent)

        axs[end].xticks = niceticks(vcat(x1mth, x2mth), 8)
        axs[end].yticks = niceticks(vcat(y1mth, y2mth), 8)

        row, col = nextpos(row, col, maxcols)
    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)
    fig
end


function exists_met(est)
    colnames = Ref(names(est))
    fields = ["tmin", "tmax", "wind", "rain", "totrad", "drrad", "dfrad"]
    lst = fields .∈ colnames
    (; tmin=lst[1], tmax=lst[2], wind=lst[3], rain=lst[4],
       totrad=lst[5], drrad=lst[6], dfrad=lst[7])
end


function plottitle(site, suffix)
    loc = replace(site, "-" => " ", "_" => " ")
    pre = isempty(site) ? "" : "$(loc): $(suffix)"
    function addtxt(txt)
        isempty(pre) ? "" : pre * txt
    end
    addtxt
end


function plot_daily(obs, est; plot::Bool=true, site::AbstractString="")
    has = exists_met(est)
    title = plottitle(site, "Daily ")
    f1 = has.tmin ? plot_daily(obs, est, :tmin, title("Tmin"), plot) : nothing
    f2 = has.tmax ? plot_daily(obs, est, :tmax, title("Tmax"), plot) : nothing
    f3 = has.wind ? plot_daily(obs, est, :wind, title("Wind"), plot) : nothing
    f1, f2, f3
end


function plot_monthly(obs, est; plot::Bool=true, site::AbstractString="")
    has = exists_met(est)
    title = plottitle(site, "Monthly ")
    f1 = has.tmin ? plot_monthly(obs, est, mean, :tmin, title("Tmin"), plot) : nothing
    f2 = has.tmax ? plot_monthly(obs, est, mean, :tmax, title("Tmax"), plot) : nothing
    f3 = has.wind ? plot_monthly(obs, est, mean, :wind, title("Wind"), plot) : nothing
    f4 = has.rain ? plot_monthly(obs, est, sum, :rain, title("Rain"), plot) : nothing
    f1, f2, f3, f4
end


function plot_kde(obs, est; plot::Bool=true, site::AbstractString="")
    has = exists_met(est)
    title = plottitle(site, "Daily ")
    f1 = has.tmin ? plot_kde(obs, est, :tmin, title("Tmin"), plot) : nothing
    f2 = has.tmax ? plot_kde(obs, est, :tmax, title("Tmax"), plot) : nothing
    f3 = has.wind ? plot_kde(obs, est, :wind, title("Wind"), plot) : nothing
    if has.rain
        # use only wet days
        obsrain = obs[obs.rain .> 0, :]
        estrain = est[est.rain .> 0, :]
        f4 = plot_kde(obsrain, estrain, :rain, title("Rain"), plot)
    else
        f4 = nothing
    end
    f1, f2, f3, f4
end


function plot_r(obs, est; plot::Bool=true, site::AbstractString="")
    has = exists_met(est)
    title = plottitle(site, "Monthly ")
    if has.tmax && has.rain
        f1 = plot_r(obs, est, [:tmax, :rain], [mean, sum], title("Tmax and Rain"), plot)
    else
        f1 = nothing
    end
    f1
end


function plot_solar(est; plot::Bool=true, site::AbstractString="")
    has = exists_met(est)
    if has.totrad && has.drrad && has.dfrad
        title1 = plottitle(site, "Daily ")
        f1 = plot_solar_daily(est, title1("Solar Radiation"), plot)
        title2 = plottitle(site, "Monthly ")
        f2 = plot_solar_monthly(est, title2("Solar Radiation"), plot)
        f3 = plot_solar_kde(est, title1("Solar Radiation"), plot)
    else
        f1 = f2 = f3 = nothing
    end
    f1, f2, f3
end


function plot_solar_daily(est, title::AbstractString, plot::Bool=true)
    gest = groupby(est, :year)

    N = length(gest)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)

    axs = []
    row = col = 1
    for i ∈ 1:N
        gi = gest[i]
        x = collect(1:length(gi.year))
        year = gi.year[1]
        y1 = gi.totrad
        y2 = gi.drrad
        y3 = gi.dfrad

        push!(axs, Axis(fig[row, col], xgridvisible=false, ygridvisible=false,
                        xlabel="DOY", ylabel="solar radiation", title=string(year)))

        s1 = scatter!(axs[end], x, y1, color=:black, marker=:circle, markersize=7)
        s2 = scatter!(axs[end], x, y2, color=:blue, marker=:xcross, markersize=7)
        s3 = scatter!(axs[end], x, y3, color=:red, marker=:utriangle, markersize=7)

        if i == N
            elem1 = MarkerElement(color = :black, marker = :circle, markersize = 15)
            elem2 = MarkerElement(color = :blue, marker = :xcross, markersize = 15)
            elem3 = MarkerElement(color = :red, marker = :utriangle, markersize = 15)
            Legend(fig[row+1, :], [elem1, elem2, elem3], ["total", "direct", "diffuse"],
                   orientation=:horizontal, framevisible=false,
                   colgap=10, backgroundcolor=:transparent)
        end

        mthdays = cummulative_days(year)
        axs[end].xticks = mthdays[1:2:end]
        axs[end].yticks = niceticks(0.0, maximum(vcat(y1, y2, y3)), 6)

        row, col = nextpos(row, col, maxcols)
    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)
    fig
end


function plot_solar_monthly(est, title::AbstractString, plot::Bool=true)
    gest = groupby(est, :year)

    N = length(gest)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)

    sbm = summarize_by_month
    axs = []
    row = col = 1
    for i ∈ 1:N
        gi = gest[i]
        x = collect(1:12)
        year = gi.year[1]
        y1 = sbm(year, Vector(gi.totrad), mean)
        y2 = sbm(year, Vector(gi.drrad), mean)
        y3 = sbm(year, Vector(gi.dfrad), mean)

        push!(axs, Axis(fig[row, col], xgridvisible=false, ygridvisible=false,
                        xlabel="Month", ylabel="solar radiation", title=string(year)))

        s1 = scatter!(axs[end], x, y1, color=:black, marker=:circle, markersize=12)
        s2 = scatter!(axs[end], x, y2, color=:blue, marker=:xcross, markersize=12)
        s3 = scatter!(axs[end], x, y3, color=:red, marker=:utriangle, markersize=12)

        if i == N
            elem1 = MarkerElement(color = :black, marker = :circle, markersize = 15)
            elem2 = MarkerElement(color = :blue, marker = :xcross, markersize = 15)
            elem3 = MarkerElement(color = :red, marker = :utriangle, markersize = 15)
            Legend(fig[row+1, :], [elem1, elem2, elem3], ["total", "direct", "diffuse"],
                   orientation=:horizontal, framevisible=false,
                   colgap=10, backgroundcolor=:transparent)
        end

        axs[end].xticks = [1:2:12...]
        axs[end].yticks = niceticks(0.0, maximum(vcat(y1, y2, y3)), 6)

        row, col = nextpos(row, col, maxcols)
    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)
    fig
end


function plot_solar_kde(est, title::AbstractString, plot::Bool=true)
    gest = groupby(est, :year)

    N = length(gest)
    maxcols = 3
    fig, nrows, ncols = create_fig(N, maxcols, title)

    axs = []
    row = col = 1
    for i ∈ 1:N
        gi = gest[i]
        x = collect(1:length(gi.year))
        year = gi.year[1]
        y1 = gi.totrad
        y2 = gi.drrad
        y3 = gi.dfrad

        push!(axs, Axis(fig[row, col], xgridvisible=false, ygridvisible=false,
                        ylabel="solar radiation", title=string(year)))

        x1 = repeat([1], length(x))
        lmt = extrema(y1)
        violin!(axs[end], x1, y1; color=:gray85, strokewidth=1, width=0.5,
                orientation=:vertical, datalimits=lmt)

        x2 = repeat([2], length(x))
        lmt = extrema(y2)
        violin!(axs[end], x2, y2; color=:gray85, strokewidth=1, width=0.5,
                orientation=:vertical, datalimits=lmt)

        x3 = repeat([3], length(x))
        lmt = extrema(y3)
        violin!(axs[end], x3, y3; color=:gray85, strokewidth=1, width=0.5,
                orientation=:vertical, datalimits=lmt)

        boxplot!(axs[end], x1, y1; orientation=:vertical, width=0.1,
                color = (:gray, 0.3), strokewidth=1, markersize=7, outliercolor=:black)

        boxplot!(axs[end], x2, y2; orientation=:vertical, width=0.1,
                 color=(:gray, 0.3), strokewidth=1, markersize=7, outliercolor=:black)

        boxplot!(axs[end], x3, y3; orientation=:vertical, width=0.1,
                 color=(:gray, 0.3), strokewidth=1, markersize=7, outliercolor=:black)

        axs[end].xticks = (1:3, ["total", "direct", "diffuse"])
        row, col = nextpos(row, col, maxcols)

    end

    foreach(c -> colsize!(fig.layout, c, Aspect(1, 1.0)), 1:ncols)
    resize_to_layout!(fig)
    plot && display(fig)
    fig
end
