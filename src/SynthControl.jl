module SynthControl

using Dates, DataFrames, Optim, LinearAlgebra, RecipesBase, Parameters
import CSV

using Reexport
@reexport using TreatmentPanels

export SynthControlModel, fit!, isfitted, load_brexit, load_germany, load_basque, load_smoking

include("synthcontrolmodels.jl")
include("plotrecipes.jl")
include("load_data.jl")

end # module
