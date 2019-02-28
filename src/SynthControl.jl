module SynthControl

import Base:show
using DataFrames, Optim, LinearAlgebra
export estimate_model

df = CSV.read("C:/Users/Nils-Holger/Desktop/Input data - cost of Brexit 2018q3.csv")

# Estimation target
goal = df[(df.country .== "United Kingdom") .& (df.dateid .< 86), :realgdp]

# Get dimensions of problem
no_comps = length(unique(df.country))-1
no_pretreat_t = length(unique(df[df.dateid .< 86, :].dateid))

# Create comparator matrix
comps = df[.!(df.country .== "United Kingdom") .& (df.dateid .< 86), :realgdp]
comps = collect(reshape(comps, no_comps, no_pretreat_t)')

# Objective function
obj(w) = norm(goal - sum((w .* comps')', dims = 2)) + 1_000_000*abs(1-sum(w))

# Initial condition and bounds
initial = rand(no_comps)
lower = zeros(no_comps)
upper = ones(no_comps)

# Get weights through optimization
res = optimize(obj, lower, upper, zeros(no_comps))
w_res = res.minimizer

# Plot result
using Plots
plot(df[(df.country .== "United Kingdom"), :realgdp], label = "UK", legend = :topleft)
plot!(sum((w_res .* reshape(df[.!(df.country .== "United Kingdom"), :realgdp],
                no_comps, length(unique(df.dateid))))', dims = 2), label = "Control")
vline!([no_pretreat_t], label = "Brexit")

end # module
