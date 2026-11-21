
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


# A constant series: its spread is no more than rounding noise (a root-mean-square
# deviation under 1e-9 of the mean's size). Float rounding leaves a constant series like
# fill(0.3, 31) a tiny spread, enough to give it a skew of ±1 and a lag-1 correlation of
# 0.97 unless it is caught first.
const FLAT_TOL = 1e-9
isflat(s2, n, m) = sqrt(s2 / n) <= FLAT_TOL * max(abs(m), 1.0)    # s2: sum of squared deviations
isflat(x::AbstractVector) = (m = mean(x); isflat(sum(v -> (v - m)^2, x), length(x), m))


function acf1(data::AbstractVector)
    isflat(data) && return 0.0
    r = first(autocor(data, [1]))
    isnan(r) ? 0.0 : r
end


# Skewness, taken as 0 for a constant series (as the weather designer does)
function skew0(data::AbstractVector)
    isflat(data) && return 0.0
    s = skewness(data)
    isnan(s) ? 0.0 : s
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
    flat = isflat(s2, n, m)
    (mean=m, sd=n > 1 ? sqrt(s2 / (n - 1)) : 0.0, rlag=flat ? 0.0 : lag / s2,
     skew=flat ? 0.0 : (s3 / n) / sqrt((s2 / n)^3))
end


# The latitude on a file's first line, or nothing if it is not a number. Excel may save
# the line with trailing commas ("3.0,,,,") and a byte-order mark first; both are dropped.
parse_lat(line::AbstractString) =
    tryparse(Float64, strip(replace(strip(line), r"^\ufeff" => "", r"[,;\s]+$" => "")))


function csv2df(fname::AbstractString)
    open(fname, "r") do fin
        line = readline(fin)
        lat = parse_lat(line)
        isnothing(lat) && error("the first line of $(basename(fname)) must be the latitude " *
                                "(a number), not \"$(first(line, 40))\"")
        df = DataFrame(CSV.File(fin; comment="#", ignoreemptyrows=true))
        (; df, lat)
    end
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
