using Test
using SynthControl, CSV, DataFrames, Dates

@testset "brexit" begin
  df = load_brexit()

  tp = TreatmentPanel("United Kingdom" => Date(2016, 7, 1), df;
                       outcome = :realgdp, id_var = :country, t_var = :quarter)

  s = SynthControlModel(tp)

  @test !(isfitted(s))

  fit!(s)

  @test isfitted(s)

  @test s.tp.no_comps == 22

  @test s.tp.no_pretreat_t == 29

  @test size(s.tp.comps) == (22, 30)

end
