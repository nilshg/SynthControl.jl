using Test
using SynthControl, CSV, Plots, DataFrames

@testset "brexit" begin
  df = load_brexit()

  tp = TreatmentPanel(data = df, outcome = :realgdp, id_var = :country, t_var = :quarter,
                      treatment = "United Kingdom, 86")

  s = SynthControlModel(tp)

  @test !(isfitted(s))

  fit!(s)

  @test isfitted(s)

  @test s.tp.no_comps == 22

  @test s.tp.no_pretreat_t == 29

  @test size(s.tp.comps) == (22, 30)

end
