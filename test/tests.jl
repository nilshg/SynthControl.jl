using Revise
using SynthControl, CSV, Plots, DataFrames

df = load_brexit()

s = SynthControlModel(df, :country, :dateid, :realgdp, 86, "United Kingdom")

@assert isfitted(s) == false

fit!(s)

@assert isfitted(s) == true

@assert s.no_comps == 22

plot(s)
