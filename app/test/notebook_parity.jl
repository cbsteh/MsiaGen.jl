# Checks that src/designer.jl gives exactly what the Pluto notebook ui/weather_designer.jl
# gives, for 40 random designs (stats file, illustrative days, rain strips) and for the
# average year of an observed file. Needs the notebook packages (run in a scratch env):
#   julia --project=<env with Distributions, HypertextLiteral, PlutoUI, SpecialFunctions> test/notebook_parity.jl
# Retire with the notebook.

using Random
include(joinpath(@__DIR__, "..", "..", "ui", "weather_designer.jl"))   # notebook as a script, defaults
module D
using Dates, Distributions, Random, SpecialFunctions, Statistics
include(joinpath(@__DIR__, "..", "src", "designer.jl"))
end
rng = MersenneTwister(7)
pick(r) = r[rand(rng, 1:length(r))]
nfail = 0
for trial in 1:40
    # random settings, built for both
    tw = Dict{Symbol,Any}(); tt = Dict{Symbol,Any}(); tm = Dict{Symbol,Any}()
    dz = Dict{String,Any}()
    for v in TVARS
        whole = Dict(k => pick(getproperty(T_RANGE, k)) for k in (:M,:A,:peak,:sd,:rlag,:skew))
        targ = Dict(k => pick(getproperty(T_RANGE, k)) for k in (:M,:A,:peak,:sd,:rlag,:skew))
        path = pick(["linear","accel","level","scurve"])
        mon = Dict(k => [pick(getproperty(T_RANGE, r)) for _ in 1:12] for (k, r) in (("dm",:dm),("sd",:dsd),("rl",:drl),("sk",:dsk)))
        for (k,x) in whole; tw[Symbol("$(v)_$k")] = x; end
        for (k,x) in targ; tt[Symbol("$(v)_$k")] = x; end
        tt[Symbol("$(v)_path")] = path
        for (k,xs) in mon, m in 1:12; tm[Symbol("$(v)_$k$m")] = xs[m]; end
        dz[v] = Dict{String,Any}("whole" => Dict(String(k)=>x for (k,x) in whole),
            "target" => merge(Dict{String,Any}("path"=>path), Dict(String(k)=>targ[k]-whole[k] for k in keys(targ))),
            "month" => Dict("dm"=>mon["dm"],"dsd"=>mon["sd"],"drl"=>mon["rl"],"dsk"=>mon["sk"]))
    end
    ww = (; (k => pick(getproperty(W_RANGE, k)) for k in (:M,:A,:peak,:sd,:rlag))...)
    wt = (; path=pick(["linear","accel","level","scurve"]), (k => pick(getproperty(W_RANGE, k)) for k in (:M,:A,:peak,:sd,:rlag))...)
    wmon = Dict(k => [pick(getproperty(W_RANGE, r)) for _ in 1:12] for (k, r) in (("dm",:dm),("sd",:dsd),("rl",:drl)))
    wm = (; (Symbol("$k$m") => wmon[k][m] for k in keys(wmon) for m in 1:12)...)
    dz["wind"] = Dict{String,Any}("whole"=>Dict(String(k)=>getfield(ww,k) for k in keys(ww)),
        "target"=>merge(Dict{String,Any}("path"=>wt.path), Dict(String(k)=>getfield(wt,k)-getfield(ww,k) for k in keys(ww))),
        "month"=>Dict("dm"=>wmon["dm"],"dsd"=>wmon["sd"],"drl"=>wmon["rl"]))
    rw = (; (k => pick(getproperty(R_RANGE, k)) for k in (:total,:seas,:dry,:pww,:pwd))...)
    rt = (; path=pick(["linear","accel","level","scurve"]), (k => pick(getproperty(R_RANGE, k)) for k in (:total,:seas,:dry,:pww,:pwd))...)
    rmon = Dict(k => [pick(getproperty(R_RANGE, r)) for _ in 1:12] for (k, r) in (("dr",:dr),("ww",:dww),("wd",:dwd)))
    rm_ = (; (Symbol("$k$m") => rmon[k][m] for k in keys(rmon) for m in 1:12)...)
    dz["rain"] = Dict{String,Any}("whole"=>Dict(String(k)=>getfield(rw,k) for k in keys(rw)),
        "target"=>merge(Dict{String,Any}("path"=>rt.path),
                        Dict(String(k)=>(k == :total ? 100*(rt.total/rw.total - 1) : getfield(rt,k)-getfield(rw,k)) for k in keys(rw))),
        "month"=>Dict("dr"=>rmon["dr"],"dww"=>rmon["ww"],"dwd"=>rmon["wd"]))
    a = stats_csv(lat, years, all_t_stats(years, (; tw...), (; tt...), (; tm...)),
                  all_w_stats(years, ww, wt, wm), all_r_stats(years, rw, rt, rm_))
    st = D.all_stats(dz, years)
    b = D.stats_csv(lat, years, st)
    if a != b
        global nfail += 1
        la, lb = split(a, '\n'), split(b, '\n')
        i = findfirst(i -> la[i] != lb[i], eachindex(la))
        fa, fb = split(la[i], ','), split(lb[i], ',')
        j = findfirst(j -> fa[j] != fb[j], eachindex(fa))
        println("trial $trial differs: line $i col $(D.STATS_HEADER[j]) $(fa[j]) vs $(fb[j])")
    end
    # illustrative points and checks
    s1, s2 = all_t_stats(years, (; tw...), (; tt...), (; tm...))["tmax"][end], st.t["tmax"][end]
    t_sample_year(s1, t_draws) == D.t_sample_year(s2, D.t_draws(0)) || println("t sample differs $trial")
    w1, w2 = all_w_stats(years, ww, wt, wm)[end], st.w[end]
    w_sample_year(w1, w_draws) == D.w_sample_year(w2, D.w_draws(0)) || println("w sample differs $trial")
    r1, r2 = all_r_stats(years, rw, rt, rm_)[end], st.r[end]
    make_strip(r1, r_draws) == D.make_strip(r2, D.r_draws(0)) || println("r strip differs $trial")
end
println("csv mismatches: $nfail / 40")
# baseline from a real obs file
txt = read(joinpath(@__DIR__, "..", "..", "data", "Serdang", "Serdang-obs.csv"), String)
b1 = average_year(parse_weather(txt)); b2 = D.average_year(D.parse_weather(txt))
println("baseline equal: ", b1 == b2)
d = D.temp_design("tmax", b2.tmax)
println("T init: ", d["whole"], " ", d["month"]["dm"])
