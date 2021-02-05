"""
    TreatmentPanel

Holds the data to be used in fitting synthetic control models. The object is constructed
from a `DataFrame` which contains all outcome and covariate data. The constructor requires
passing information on which column in the `DataFrame` holds the subject and time period
identifiers, as well as information on the timing of treatment. Treatments are specified
as a `Pair` of either a `String` or `Symbol` identifying the unit treated, and a `Date` or
`Int` value, identifying the treatment period. Where treatments have an end point, the
treatment period os specified as a `Tuple` of either `Date` or `Int`, indicating the first
and last period of treatment. Where individual units have multiple treatment periods,
these are specified as a `Pair{Union{String, Symbol}, Vector{Union{Int, Date}}}` of treatment
unit identifier and vector of timings. Finally, single continuous, single periodic, and
multiple period treatments can be generalized to multiple treated units.

The following table provides an overview of the types of treatment pattern supported:

|                 |  Only starting point        |   Start and end point                     |   Multiple start & end points                       |
|-----------      |------------------------     |-------------------------                  |------------------------------                       |
| **one unit**         |  Pair{String, Date}         |   Pair{String, Tuple{Date, Date}}         |  Pair{String}, Vector{Tuple{Date, Date}}}           |
| **multiple units**  |  Vector{Pair{String, Date}} |   Vector{Pair{String, Tuple{Date, Date}}} |  Vector{Pair{String}, Vector{Tuple{Date, Date}}}}   |

Currently, only single treatment unit and continuous treatment is supported.

"""

@with_kw struct TreatmentPanel
    data::DataFrame
    outcome::Union{Symbol, String}
    id_var::Union{Symbol, String}
    t_var::Union{Symbol, String}
    treatment::Union{Pair{Union{Symbol, String}, Union{Date, Int64}}}
    y₁₀::Vector{Float64} # (T₀ × 1) vector of pretreatment outcomes for treated unit
    y₁₁::Vector{Float64} # (T-T₀ × 1) vector of posttreatment outcomes for treated unit
    yⱼ₀::Matrix{Float64} # (T₀ × J) matrix of pretreatment outcomes for donor pool
    yⱼ₁::Matrix{Float64} # (T-T₀ × J) matrix of posttreatment outcomes for donor pool
    x₁₀::Union{Nothing, Vector{Float64}} # (k × 1) vector of predictors for treated unit
    xⱼ₀::Union{Nothing, Matrix{Float64}} # (k × J) matrix of predictors for donor pool
    predictors::Union{Nothing, Union{Symbol, String}, Vector{Union{Symbol, String}}}
    J::Int64
    T₁::Int64
    T₀::Int64
    comp_labels::Vector{Union{Symbol, String}}
end

function check_tp(outcome, id_var, t_var, treatment, data)

    treatᵢ, treatₜ = treatment

    !isnothing(outcome) || error(ArgumentError(
        "Please specify outcome, the nameof the column in your dataset holding the "*
        "outcome variable of interest in your dataset."
        ))
    !isnothing(id_var) || error(ArgumentError(
            "Please specify id_var, the name of the column in your data set holding the "*
            "identifier of your units of observation (panel dimension)."
            ))
    !isnothing(t_var) || error(ArgumentError(
            "Please specify t_var, the name of the column in your dataset holding the "*
            "time dimension."
        ))
    !isnothing(treatment) || error(ArgumentError(
        "Please specify treatment, the identifier of the treated unit(s) in your "*
        "dataset and associated start(s) of treatment."
        ))
    in(string(outcome), names(data)) || error(DomainError(
        "The specified outcome column, $outcome, was not found in your data set. Please "*
        "double check the name of the outcome column."
        ))
    in(treatᵢ, data[!, id_var]) || error(DomainError(
        "The specified treatment unit, $(treatᵢ), was not found in the specified id "*
        "column, $id_var. Please ensure that treatment unit and id column are correctly "*
        "specified."
        ))
    in(treatₜ, data[!, t_var]) || error(DomainError(
        "The specified treatment period, $(treatₜ), was not found in the specified time "*
        "column, $t_var. Please double check treatment time and time column."
        ))
    treatₜ != maximum(data[!, t_var]) || error(DomainError(
        "The specified treatment period, $(treatₜ), is equal to the last observation period "*
        "of your data set - no estimation of treatment effects is possible. Please double "*
        "check the treatment period."
        ))
    treatₜ < maximum(data[!, t_var]) || error(DomainError(
        "The specified treatment period, $(treatₜ), occurs after the last observation in your "*
        "data set, $(maximum(data[!, t_var])) - no estimation of treatment effects is possible. "*
        "Please double check specification of treatment period."
        ))
    # To add - checks around balance of panel?
end

function TreatmentPanel(treatment::T, data;
                        outcome = nothing, id_var = nothing, t_var = nothing,
                        predictors = nothing) where T <: Pair

    check_tp(outcome, id_var, t_var, treatment, data)

    treatᵢ, treatₜ = treatment

    y₁₀ = data[(data[!, id_var] .== treatᵢ) .&
               (data[!, t_var] .< treatₜ), outcome]
    yⱼ₀ = data[(data[!, id_var] .!= treatᵢ) .&
               (data[!, t_var] .< treatₜ), outcome]

    y₁₁ = data[(data[!, id_var] .== treatᵢ) .&
               (data[!, t_var] .>= treatₜ), outcome]
    yⱼ₁ = data[(data[!, id_var] .!= treatᵢ) .&
               (data[!, t_var] .>= treatₜ), outcome]

    T₁ = length(unique(data[data[!, t_var] .>= treatₜ, t_var]))
    T₀ = length(y₁₀)
    J = length(unique(data[data[!, id_var] .!= treatᵢ, id_var]))

    comp_labels = unique(data[data[!, id_var] .!= treatᵢ, id_var])

    if predictors == nothing
        TreatmentPanel(
            data,
            outcome,
            id_var,
            t_var,
            treatment,
            y₁₀,
            y₁₁,
            reshape(yⱼ₀, J, T₀),
            reshape(yⱼ₁, J, T₁),
            nothing,
            nothing,
            nothing,
            J,
            T₁,
            T₀,
            comp_labels
        )
    else
        "Not yet implemented"
    end
end
