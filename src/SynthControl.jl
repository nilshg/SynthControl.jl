module SynthControl

using Dates, DataFrames, LinearAlgebra, RecipesBase, Statistics
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

# Utilities
include("plotrecipes.jl")
include("load_data.jl")

# Precompilation
using PrecompileTools

@setup_workload begin
    df = DataFrame(state = [fill("California", 5); fill("Wyoming", 5)],
            year = [1980:1984; 1980:1984],
            outcome = [rand(4); rand()+3; rand(5)])

    @compile_workload begin
        bp = BalancedPanel(df, 
            "California" => 1984; id_var = :state, t_var = :year, outcome_var = :outcome)

        simplescm = SimpleSCM(bp)
        fit!(simplescm)

        synthdid = SyntheticDiD(bp)
        fit!(synthdid)
    end
end

end # module
