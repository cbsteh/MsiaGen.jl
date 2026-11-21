# Start the MsiaGen UI:  julia --project=. -t auto start.jl [port]
# (or double-click StartMsiaGen.command on macOS, StartMsiaGen.bat on Windows)
using Sockets
port = isempty(ARGS) ? 8100 : parse(Int, ARGS[1])

# stop at once (before the slow package loading) if MsiaGen already runs in another window
try
    close(listen(Sockets.localhost, port))
catch
    println("\nMsiaGen is already running in another window (port $port).",
            "\nClose that window first, then start MsiaGen again.\n")
    exit(1)
end

# Ctrl+C stops the app quietly: by default a script dies on Ctrl+C with a dump of every
# thread's stack; this turns it into an exception caught below
Base.exit_on_sigint(false)
try
    println("  Loading MsiaGen...")
    include(joinpath(@__DIR__, "src", "WeatherApp.jl"))
    # (invokelatest: the module was defined inside this same block)
    Base.invokelatest(getfield(Main, :WeatherApp).start; port, open_browser = get(ENV, "MSIAGEN_NO_BROWSER", "") == "")
    while true
        sleep(1)
    end
catch e
    e isa InterruptException || rethrow()
    println("\n  MsiaGen stopped.")
    exit(0)
end
