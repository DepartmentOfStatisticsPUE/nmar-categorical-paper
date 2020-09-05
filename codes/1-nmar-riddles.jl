## riddles et al nmar

using Distributions
using Random
using DataFrames
using StatsBase
using Statistics
using FreqTables

strata_names = [11,12,21,22]
strata_sizes = 500 .* [50,100,150,200]
## x = vcat(fill.(strata_names, strata_sizes)...)
x = inverse_rle(strata_names, strata_sizes)
