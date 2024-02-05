module SynthControl

using Dates, DataFrames, LinearAlgebra, RecipesBase, Statistics, StatsBase
using JuMP, HiGHS 

import CSV

# Re-export TreatmentPanels for easier interaction with the input data
using Reexport
@reexport using TreatmentPanels

# Export main functionality
export SimpleSCM, SyntheticDiD, fit!, isfitted

# Export example data utilities
export load_brexit, load_germany, load_basque, load_smoking
export load_brexit_panel, load_germany_panel, load_basque_panel, load_smoking_panel

# Estimators
include("SimpleSCM.jl")
include("SyntheticDiD.jl")
include("MC_fect.jl")

# Utilities
include("plotrecipes.jl")
include("load_data.jl")

# Precompilation
using PrecompileTools

@setup_workload begin
    @compile_workload begin
        sd = load_smoking()
        sp = load_smoking_panel()

        simplescm = SimpleSCM(sp)
        fit!(simplescm)

        synthdid = SyntheticDiD(sp)
        fit!(synthdid)

        fect_result = SynthControl.fect_default(sp; verbose = false)
    end
end

end # module
