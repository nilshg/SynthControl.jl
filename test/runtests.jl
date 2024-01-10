using Test
using SynthControl, TreatmentPanels, DataFrames, Dates

@testset "SimpleSCM" begin

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
  @test tp isa BalancedPanel{SingleUnitTreatment{Continuous}}

  s = SimpleSCM(tp)

  @test !(isfitted(s))

  fit!(s)

  @test isfitted(s)

end

@testset "Synthetic Differences-in-Differences" begin

  bp = load_smoking_panel()

  @test bp isa BalancedPanel{SingleUnitTreatment{Continuous}}

  sdid_model = SyntheticDiD(bp)

  fit!(sdid_model)

  @test only(sdid_model.τ̂) ≈ -15.604 atol=0.1

  @test_throws "placebo" fit!(sdid_model, se = :jackknife)

  fit!(sdid_model, se = :placebo)

  @test only(sdid_model.se_τ̂) ≈ 9.3 atol=0.1
end


@testset "Data sets" begin

  df = load_germany()
  @test only(df[df.country .== "USA" .&& df.year .== 1961, :gdp]) == 2929

  bp = load_germany_panel()
  @test typeof(bp) <: BalancedPanel

  df = load_basque()
  @test only(df[df.regionno .== 1 .&& df.year .== 1957, :gdpcap]) == 2.60361314875437

  bp = load_basque_panel()
  @test typeof(bp) <: BalancedPanel

  df = load_brexit()
  @test only(df[df.country .== "Australia" .&& df.quarter .== Date(2009, 1, 1), :realgdp]) == 1.04

  bp = load_brexit_panel()
  @test typeof(bp) <: BalancedPanel

end