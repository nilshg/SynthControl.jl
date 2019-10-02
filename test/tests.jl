using Test
using SynthControl, CSV, Plots, DataFrames

@testset "brexit" begin
  df = load_brexit()

  s = SynthControlModel(df, :country, :dateid, :realgdp, 86, "United Kingdom")

  @assert isfitted(s) == false

  fit!(s)

  @assert isfitted(s) == true

  @assert s.no_comps == 22
end


