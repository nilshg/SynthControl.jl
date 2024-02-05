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

@testset "SyntheticDiD" begin

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


@testset "MC-NNM" begin
  # Test outputs of helper functions against R/C++ versions
  test22 = [1.0 3; 2 4]
  test32 = [1.0 4; 2 5; 3 6]
  test23 = [1.0; 2;; 3; 4;; 5; 6]

  YYt, μ_Yt, α_Yt, xᵢ_Yt = SynthControl.Y_demean(test22, 1)
  @test all((YYt, μ_Yt, α_Yt, xᵢ_Yt) .== ([-0.5 -0.5; 0.5 0.5], 2.5, [1.5, 3.5], [0.0, 0.0]))
  YYt, μ_Yt, α_Yt, xᵢ_Yt = SynthControl.Y_demean(test22, 2)
  @test all((YYt, μ_Yt, α_Yt, xᵢ_Yt) .== ([-1.0 1.0; -1.0 1.0], 2.5, [0.0, 0.0], [2.0, 3.0]))
  YYt, μ_Yt, α_Yt, xᵢ_Yt = SynthControl.Y_demean(test22, 3)
  @test all((YYt, μ_Yt, α_Yt, xᵢ_Yt) .== ([0 0; 0 0], 2.5, [1.5, 3.5], [2.0, 3.0]))

  μ_test, FE_ad_test, α_test, xᵢ_test = SynthControl.fe_add([1, 2], [3, 4], 5, 2, 2, 3)
  @test all((μ_test, FE_ad_test, α_test, xᵢ_test) .== (5, [-1.0 0.0; 0.0 1.0], [-4, -3], [-2, -1]))

  μt, FEt, αt, xᵢt, FEiut = SynthControl.ife(test22, 3, 1, 1, false, 0.1)
  @test all(isapprox.((μt, FEt, αt, xᵢt, FEiut), (2.5, [1.0 3; 2 4], [-1.0, 1], [-0.5, 0.5], [0.0 0; 0 0])))
  μt, FEt, αt, xᵢt, FEiut = SynthControl.ife(test32, 3, 1, 1, false, 0.1)
  @test all(isapprox.((μt, FEt, αt, xᵢt, FEiut), (3.5, [1.0 4; 2 5; 3 6], [-1.5, 1.5], [-1.0, 0, 1], [0.0 0; 0 0; 0 0])))

  λt, ft, VNTt, FEt = SynthControl.panel_factor(test22, 1)
  @test all(isapprox.((λt, ft, VNTt, FEt), ([-0.5721252, -1.2933185], [-2.226040, -3.158762], [7.466517], 
      [1.273574 2.878979; 1.807207 4.085286]), atol = 0.0001))
  λt, ft, VNTt, FEt = SynthControl.panel_factor(test22, 2)
  @test all(isapprox.((λt, ft, VNTt, FEt), ([-0.5721252 -1.2933185; -1.2933185 0.5721252], 
      [-2.226040 0.2115285; -3.158762 -0.1490682], [7.466517 0; 0 0.03348281], test22), atol = 0.01))

  @test SynthControl.panel_FE(test22, 0.1, false) ≈ [1.180357 2.668257; 1.674932 3.786270] atol = 0.001
  @test SynthControl.panel_FE(test32, 0.1, false) ≈ [1.346675 3.575952; 1.930930 4.660546; 2.515185 5.745139] atol = 0.001
  @test SynthControl.panel_FE(test23, 0.1, false) ≈ [1.271176 2.902109 4.533041; 1.610219 3.676147 5.742075] atol = 0.001
  @test SynthControl.panel_FE(test23, 0.2, false) ≈ [1.185724 2.707020 4.228317; 1.501975 3.429025 5.356076] atol = 0.001
  #!# hard impute = false not tested as the function does not seem exported in the R package?

  @test SynthControl.fect_default(load_smoking_panel(), verbose = false).τ̂ᵃᵗᵗ ≈ -25.907 atol = 0.01
end