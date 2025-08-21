
macro m(idx)
    quote
        $(esc(idx)) + 1
    end
end


function days_in_each_month(year::Int)
    feb =  isleapyear(year) ? 29 : 28
    [31, feb, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
end


cummulative_days(year::Int) = accumulate(+, days_in_each_month(year))


function acf1(data::AbstractVector)
    r = first(autocor(data, [1]))
    isnan(r) ? 0 : r
end


first_row_to_vec(df::DataFrame) = Vector(df[1, :])


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
    f_lst = collect(Iterators.flatten(lst))
    txt = "%8.2f " ^ length(f_lst)
    txt = prefix * txt * "\n"
    fmt = Printf.Format(txt)
    Printf.format(stdout, fmt, f_lst...)
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


function check_errors(thd::AbstractVector, tgt::AbstractVector, est::AbstractVector;
                      p=0.99)
    f_thd = collect(Iterators.flatten(thd))
    f_tgt = collect(Iterators.flatten(tgt))
    f_est = collect(Iterators.flatten(est))
    error = 100 * abs.(f_est - f_tgt) ./ abs.(f_tgt)
    delta = 100 * (error - f_thd) ./ f_thd
    p_ok = count(delta .<= 0) / length(delta)
    allok = p_ok >= p
    allok, error, delta
end


function summarize_by_month(year::Int, ar::Vector{Float64}, fn::Function)
    monthdays = days_in_each_month(year)
    mthx = zeros(Float64, 12)
    t0 = 1
    # iterate each month and collect its corresponding summarized data:
    for (i, num) ∈ enumerate(monthdays)
        t1 = t0 + num - 1      # work on data from the first to last day of every month
        arx = ar[t0:t1]        # slice the data according to the month
        mthx[i] = fn(arx)      # call the summary function and its arguments
        t0 = t1 + 1            # first day of next month
    end
    mthx
end


function mean_error(delta)
    delta0 = delta[delta .> 0]
    isempty(delta0) ? mean(delta) : mean(delta0)
end


function partition_list(list::AbstractVector, lengths::AbstractVector)
    partitions = Vector{Vector}()  # an empty vector of vectors
    start_index = 1
    for len ∈ lengths
        end_index = start_index + len - 1
        push!(partitions, list[start_index:end_index])
        start_index = end_index + 1
    end
    partitions
end
