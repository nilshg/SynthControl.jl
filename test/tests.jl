using Revise
using SynthControl, CSV, Plots, DataFrames

df = CSV.read("C:/Users/Nils-Holger/Desktop/Input data - cost of Brexit 2018q3.csv")
df = dropmissing(df[:, 1:4], disallowmissing=true)

s = SynthControlModel(df, :country, :dateid, :realgdp, 86, "United Kingdom")

@assert isfitted(s) == false

fit!(s)

@assert isfitted(s) == true

@assert s.no_comps == 22

plot(s)
