
macro m(idx)
    quote
        $(esc(idx)) + 1
    end
end


const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


function days_in_each_month(year::Int)
    feb =  isleapyear(year) ? 29 : 28
    [31, feb, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
end


# Day ranges of the 12 months of `year`, e.g. 1:31, 32:59, ...
function month_ranges(year::Int)
    monthdays = days_in_each_month(year)
    t1 = cumsum(monthdays)
    [(t - n + 1):t for (t, n) ∈ zip(t1, monthdays)]
end


# Views of a year's daily `data`: the whole year, then each month
function year_and_months(data::AbstractVector, year::Int)
    [view(data, r) for r ∈ [[1:length(data)]; month_ranges(year)]]
end


function acf1(data::AbstractVector)
    r = first(autocor(data, [1]))
    isnan(r) ? 0.0 : r
end


# Mean, sd, lag-1 autocorrelation and skewness of `x` (as mean, std, acf1
# and skewness give them), in two passes and without allocating
function month_stats(x::AbstractVector)
    n = length(x)
    m = mean(x)
    s2 = s3 = lag = 0.0
    d0 = 0.0
    for (k, v) ∈ enumerate(x)
        d = v - m
        d2 = d * d
        s2 += d2
        s3 += d2 * d
        k > 1 && (lag += d0 * d)
        d0 = d
    end
    rlag = lag / s2
    (mean=m, sd=sqrt(s2 / (n - 1)), rlag=isnan(rlag) ? 0.0 : rlag,
     skew=(s3 / n) / sqrt((s2 / n)^3))
end


function csv2df(fname::AbstractString)
    nt =(;)
    open(fname, "r") do fin
        lat = parse(Float64, readline(fin))
        df = DataFrame(CSV.File(fin; comment="#", ignoreemptyrows=true))
        nt = (; df, lat)
    end
    nt
end


function pprintf(lst, prefix)
    txt = prefix * "%8.2f " ^ length(lst) * "\n"
    Printf.format(stdout, Printf.Format(txt), lst...)
end


function print_start(year::Int, tgt::AbstractVector, thd::AbstractVector)
    println("Year: $year")
    pprintf(tgt, "TGT: ")
    pprintf(thd, "THD: ")
end


function print_update(ok::Bool, err)
    pprintf(err, "ERR: ")
    errtxt = ok ? "** success **" : "~ above threshold ~"
    println(errtxt)
end


# Errors (%) of the estimates `est` from the targets `tgt`, and whether at
# least a share `p` of them are within their thresholds `thd` (%). A
# target of 0 is taken as 0.001, so the error stays finite.
function check_errors(thd::AbstractVector, tgt::AbstractVector, est::AbstractVector;
                      p=0.99)
    error = 100 .* abs.(est .- tgt) ./ max.(abs.(tgt), 1e-3)
    allok = count(error .<= thd) / length(error) >= p
    allok, error
end
