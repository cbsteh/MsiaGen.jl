# MsiaGen app

A graphical interface for MsiaGen, written in Julia with Stipple, StippleUI and StipplePlotly
(Genie Framework), in the style of the Sawit app. It designs a site's weather statistics,
runs the generator and checks its fit, on the same `data/` folder that `src/main.jl` uses.

## Start it

Double-click `StartMsiaGen.command` (macOS) or `StartMsiaGen.bat` (Windows), or:

```
cd MsiaGen/app
julia --project=. -t auto start.jl [port]   # http://localhost:8100 opens in the browser
```

The first start takes a minute or two while Julia compiles. Keep the terminal window open.
The page needs no internet (fonts are in `ui/public/fonts`).

| Variable | Use |
|---|---|
| `MSIAGEN_DATA` | the data folder (default `MsiaGen/data`) |
| `MSIAGEN_NO_BROWSER` | set to anything to start without opening the browser |

## Pages

- **Sites**: the sites in the data folder, with the years of their observed, statistics and
  simulated files. **New site** designs a site from scratch.
- **1 Site**: name, latitude, start and end year, and the weather the design starts from
  (typical Malaysian values, the site's observed or simulated weather, or an uploaded file).
- **2 Temperature, 3 Wind, 4 Rainfall**: each in three tabs, *Start year* (whole-year
  values), *Each month* (changes per month) and *End year* (change by the end year and the
  path to it). One chart shows the start year and the end year together (tick boxes show or
  hide each year, the illustrative days drawn with MsiaGen's equations, and on the
  temperature page the other temperature); the table shows the monthly values and checks.
- **5 Review**: checks, annual values, the path of the design over the years, and the
  stats file: **Save to data folder** or **Download**.
- **Generate**: runs `generate_weather` from the design, from the site's stats file, or
  from its observed weather; shows MsiaGen's messages and the simulated weather by year.
- **Fit check**: `check_fit`'s comparison as a table and one 1:1 chart per statistic, with
  the fit tolerance; **Write report files** saves `-fit.txt`, `-fit.png` and `-sim.png`.

End-year targets are changes from the start-year values, so editing the start year keeps
the scenario.

## Files

```
start.jl, StartMsiaGen.command, StartMsiaGen.bat   launchers
src/WeatherApp.jl      reactive state, handlers, routes, page setup
src/designer.jl        the designer's science (seasonal curves, paths, illustrative days, checks)
src/charts.jl          Plotly charts
src/sites.jl           data folder, generating, fit check (calls MsiaGen)
ui/app.jl.html         page layout (Quasar components, Vue syntax)
ui/public/msiagen.css  styling (Sawit's colours and fonts)
```

## Limits

- Single user, local only. The design lives in the running app: it is saved only as the
  stats file. Opening a site starts its design from its observed weather (or typical values),
  not from its saved stats file.
