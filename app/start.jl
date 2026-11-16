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

println("  Loading MsiaGen...")
include(joinpath(@__DIR__, "src", "WeatherApp.jl"))
WeatherApp.start(; port, open_browser = get(ENV, "MSIAGEN_NO_BROWSER", "") == "")
while true
    sleep(1)
end
