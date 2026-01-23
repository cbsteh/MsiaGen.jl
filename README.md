
# MsiaGen: Stochastic Daily Weather Generator for Malaysia (Julia)

## Overview
MsiaGen is a stochastic daily weather generator for Malaysia’s tropical climate, emphasizing computational simplicity, site-specific parameterization, and practical applicability.

The model was calibrated using data from 12 sites across Malaysia and validated at 11 independent sites, encompassing diverse climatic conditions from Peninsular to East Malaysia. MsiaGen uses a Skew Normal distribution for air temperatures to capture observed asymmetries, particularly in maximum temperatures, while utilizing Weibull and Gamma distributions for wind speed and rainfall, respectively. The generator incorporates first-order autoregressive processes for temporal dependencies and a two-state Markov chain for wet/dry day sequencing.

Validation showed strong monthly-scale performance, with mean absolute errors below 1.2% for temperatures, 2.4% for wind speed, and 1.8% for rainfall, along with near-zero model bias and high overall model agreement scores (Kling-Gupta Efficiency metric >0.8). Daily scale validation using quantile-quantile plots revealed excellent agreement for temperature distributions, with points clustering tightly along the identity line within common ranges (21–28 °C for minimum and 25–39 °C for maximum temperatures). Empirical cumulative distribution function analysis indicated that 85±10% of daily temperature errors were within ±2.0°C, 94±6% of wind speed errors were within ±1.0 m s⁻¹, and 83±5% of rainfall errors were within ±20 mm. However, performance declined for extreme events, particularly rainfall exceeding 80–100 mm and wind speeds above 3–4 m s-1, likely due to distribution tail limitations and short observational records (3–5 years).

Further validation using oil palm yield simulations at two independent plantation sites demonstrated that generated weather reproduced temporal dynamics across multiple planting densities. MsiaGen offers a practical and data-efficient tool for tropical agricultural research.

MsiaGen is written in Julia.

## Usage
Clone the repo by

```
https://github.com/cbsteh/MsiaGen.git
```


## Example

```
using MsiaGen
using Random

function gen_weather(csv_path, seed)
    seednum = (seed < 0) ? rand(1:typemax(Int)) : seed
    Random.seed!(seednum)

    res = csv2df(csv_path)
    nt = generate_mets(res.df; verbose=false)   # set `verbose` to true for detailed run operation and fit.
    df = collate_mets(nt)
end


seed = -1  # seed number: <0 = random runs, >0 = determinstic runs
csv_path = joinpath("data", "Serdang", "data.csv")   # read data file in ndata/Serdang/data.csv
res = gen_weather(csv_path, seed)  # generated daily weather in `res`
```

## References

Sung, C. B. S., Cheah, S. S., & Appleton, D. R. (2026). A stochastic daily weather generator for perennial crop simulations in tropical Malaysia. PLOS One. (accepted).

