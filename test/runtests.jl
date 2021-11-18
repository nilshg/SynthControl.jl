using Test
using SynthControl, TreatmentPanels, DataFrames, Dates

@testset "Basic functionality" begin

  outcome_A = collect(range(1.0; stop = 3.0, length = 10))
  outcome_B = collect(range(2.0; stop = 1.0, length = 10))
  outcome_C = 0.7 .* outcome_A .+ 0.3 .* outcome_B
  outcome_C[8:end] .+= 1.0

  df = DataFrame(id = repeat(["A", "B", "C"], inner = 10), 
    time = repeat(Date(2000):Year(1):Date(2009), 3),
    value = [outcome_A; outcome_B; outcome_C])

  tp = BalancedPanel(df, "C" => Date(2007);
    outcome_var = :value, id_var = :id, t_var = :time)

  # Test that BalancedPanel works
  @test tp isa BalancedPanel{SingleUnitTreatment, ContinuousTreatment}

  s = SynthControlModel(tp)

  @test !(isfitted(s))

  fit!(s)

  @test isfitted(s)

end
