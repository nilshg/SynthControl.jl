module SynthControl

using Dates, DataFrames, Optim, LinearAlgebra, RecipesBase, Parameters, FixedEffectModels, Statistics
import CSV

using Reexport
@reexport using TreatmentPanels

export SimpleSCM, SyntheticDiD, fit!, isfitted
export load_brexit, load_germany, load_basque, load_smoking
export load_brexit_panel, load_germany_panel, load_basque_panel, load_smoking_panel

include("SimpleSCM.jl")
include("plotrecipes.jl")
include("load_data.jl")
include("SyntheticDiD.jl")

end # module
