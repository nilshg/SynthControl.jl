using Test
using SynthControl, CSV, Plots, DataFrames

@testset "brexit" begin
  df = load_brexit()

  s = SynthControlModel(df, :realgdp, :country, :dateid, "United Kingdom", 86)

  # Test keyword constructor
  s2 = SynthControlModel(df, tid = :dateid, pid = :country, outcome = :realgdp,
        treat_id = "United Kingdom", treat_t = 86)

  @test !(isfitted(s))

  fit!(s)

  @test isfitted(s)

  @test s.no_comps == 22

  @test s.no_pretreat_t == 29

  @test size(s.comps) == (29, 22)

  @test size(s2.comps) == (29, 22)
end
