# MsiaGen: Stochastic Daily Weather Generator for Malaysia (Julia)

## Overview
MsiaGen is a stochastic daily weather generator for Malaysia’s tropical climate, emphasizing computational simplicity, site-specific parameterization, and practical applicability.

The model was calibrated using data from 12 sites across Malaysia and validated at 11 independent sites, encompassing diverse climatic conditions from Peninsular to East Malaysia. MsiaGen uses a Skew Normal distribution for air temperatures to capture observed asymmetries, particularly in maximum temperatures, while utilizing Weibull and Gamma distributions for wind speed and rainfall, respectively. The generator incorporates first-order autoregressive processes for temporal dependencies and a two-state Markov chain for wet/dry day sequencing.

Validation showed strong monthly-scale performance, with mean absolute errors below 1.2% for temperatures, 2.4% for wind speed, and 1.8% for rainfall, along with near-zero model bias and high overall model agreement scores (Kling-Gupta Efficiency metric >0.8). Daily scale validation using quantile-quantile plots revealed excellent agreement for temperature distributions, with points clustering tightly along the identity line within common ranges (21–28 °C for minimum and 25–39 °C for maximum temperatures). Empirical cumulative distribution function analysis indicated that 85±10% of daily temperature errors were within ±2.0°C, 94±6% of wind speed errors were within ±1.0 m s⁻¹, and 83±5% of rainfall errors were within ±20 mm. However, performance declined for extreme events, particularly rainfall exceeding 80–100 mm and wind speeds above 3–4 m s-1, likely due to distribution tail limitations and short observational records (3–5 years).

Further validation using oil palm yield simulations at two independent plantation sites demonstrated that generated weather reproduced temporal dynamics across multiple planting densities. MsiaGen offers a practical and data-efficient tool for tropical agricultural research.

MsiaGen is written in Julia.

### Changes since the published version
The code has been revised since the validation in the paper (see References):

- **Wet-day statistics.** `pww` and `pwd` now count only the days that have a following day in the month (the standard estimator), and the number of wet days in a month is the expected number rounded to the nearest day. Previously, a fixed term of 1/30 (about one wet day per month) offset the bias of the earlier estimator. Across 23 sites, simulated wet days per month now match the observed ones with a bias of −0.01 days (previously −0.14), and `pww` and `pwd` are reproduced about 60% more closely; rainfall totals are unchanged.
- **Start of the first year.** Temperature, like wind, starts from January's mean.
- **Faster runs**, plus charts of the simulated weather (`plot_weather`), a fit report (`check_fit`) and an interactive weather designer (`ui/weather_designer.jl`).

Statistics files made with earlier versions should be rebuilt from their observed weather (see below), as their `pww` and `pwd` use the earlier estimator. With the same seed, results differ from earlier versions.

## Installation
Requires Julia 1.11 or later.

```
git clone https://github.com/cbsteh/MsiaGen.jl.git
cd MsiaGen.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The first `using MsiaGen` precompiles the package, including a short trial run, which takes a minute or so. Later sessions start quickly.

## Usage
Open `src/main.jl`, choose the site and options at the top, and run it from the project folder:

```
julia --project=.
julia> include("src/main.jl")
```

For each run, `main.jl`:

1. generates daily weather for the site (`generate_weather`), written to `data/<site>/<site>-sim.csv`, and lists any months whose generated weather missed the fit tolerance (the best of all attempts is kept);
2. charts the simulated weather on one page (`plot_weather`), saved as `<site>-sim.png`, with a summary table per year;
3. compares the simulated weather with the site's statistics (`check_fit`), saved as `<site>-fit.txt` and `<site>-fit.png`. For each statistic, the report counts the values outside the generator's fit tolerance, and it lists any years or statistics of the stats file that the simulated weather lacks.

Options in `main.jl`:

| Option | Meaning |
|---|---|
| `site` | folder in `data/` holding the site's files |
| `seed` | `-1` for a random run; a fixed number repeats a run exactly |
| `use_stats` | `true`: generate from `<site>-stats.csv`; `false`: first rebuild it from `<site>-obs.csv` |
| `dry_day`, `dry_spell`, `dry_month`, `hot_day` | thresholds for the charts (dry day, dry spell length, dry month, hot day) |

To keep editing the source files without restarting Julia, install [Revise](https://github.com/timholy/Revise.jl) in your default environment; `main.jl` loads it when it is available.

The same steps from your own code:

```julia
using MsiaGen

df = generate_weather("Serdang"; folder="data", from_stats=true, seed=1)   # a DataFrame
plot_weather("Serdang"; folder="data", dry_day=0.5)
check_fit("Serdang"; folder="data")
```

## Data files
Each site has a folder in `data/` (this repository includes `data/Serdang` as an example):

- **`<site>-obs.csv`**: observed daily weather. The first line is the site's latitude (decimal degrees); then a header with `year,month,day` and any of `tmin`, `tmax` (°C), `wind` (m/s) and `rain` (mm); then one row per day, whole years only. Lines starting with `#` are comments. Rows may be in any order, but each day of each year must appear exactly once with a value for every variable; otherwise MsiaGen stops and names the first missing, repeated or invalid date or value. In a month with a constant value, `skew` is taken as 0 (as in the weather designer).
- **`<site>-stats.csv`**: the statistics MsiaGen generates from. The first line is the latitude; then one row per year, with, for each variable, the annual value (suffix `0`) and the 12 monthly values (suffixes `1` to `12`):
  - `tmin`, `tmax`: `mean`, `sd`, `rlag` (lag-1 autocorrelation) and `skew`, e.g. `mean_tmin0`, …, `skew_tmax12`
  - `wind`: `mean`, `sd` and `rlag`, e.g. `mean_wind1`
  - `rain`: `totrain` (total, mm), `pww` (chance of a wet day after a wet day) and `pwd` (after a dry day), e.g. `totrain0`, `pww1`

  Build it from observed weather with `use_stats = false` in `main.jl` (or `create_data_file`), or design one with the weather designer. Before generating, MsiaGen checks every monthly value and stops with a list of problems if any is outside its valid range: `sd` above 0, `rlag` between −1 and 1, wind `mean` above 0, `totrain` 0 or more, `pww` and `pwd` between 0 and 1.

## Weather designer
`ui/weather_designer.jl` is a [Pluto](https://plutojl.org) notebook for designing a site's statistics, for a start year and how they change by an end year, without observed weather (or starting from an observed file). It writes a `<site>-stats.csv` for MsiaGen. To open it:

```
julia --project=.
julia> using Pluto; Pluto.run()
```

then open `ui/weather_designer.jl` from Pluto's start page.

## Tests
```
julia --project=. -e 'using Pkg; Pkg.test()'
```

## References

Teh CBS, Cheah SS, Appleton DR (2026) A stochastic daily weather generator for perennial crop simulations in tropical Malaysia. PLoS One 21(2): e0338833. https://doi.org/10.1371/journal.pone.0338833
