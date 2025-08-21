module MsiaGen

using Base64
using CairoMakie
using CSV
using DataFrames
using Dates
using Distributions
using GLM
using HypothesisTests
using LinearAlgebra
using Parameters
using PrettyTables
using Printf
using Random
using SpecialFunctions
using Statistics
using StatsBase

include("met.jl")
include("utils.jl")
include("GoF.jl")
include("collate.jl")
include("gentemp.jl")
include("genwind.jl")
include("genrain.jl")
include("gensolar.jl")
include("input.jl")
include("template/plots.jl")
include("template/render.jl")

# GoF.jl
using .GoF
export GoF

# utils.jl
export csv2df
# collate.jl
export collate_mets, collate_stats, gof, generate_mets
# gentemp.jl
export create_temp
# genwind.jl
export create_wind
# genrain.jl
export create_rain
# gentemp.jl, genwind.jl, genrain.jl
export generate!
# plots.jl
export plot_daily, plot_monthly, plot_kde, plot_r, plot_solar_daily,
       plot_solar_monthly, plot_solar_kde, plot_solar
# render.jl
export report, batch_report_to_html, report_to_html, report_to_screen,
       fit_dist_metrics, fit_monthly_metrics
# input.jl
export create_data_file

end # module MsiaGen
