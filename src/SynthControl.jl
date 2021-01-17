module SynthControl

using Dates, DataFrames, Optim, LinearAlgebra, RecipesBase, Parameters
import CSV

export TreatmentPanel, SynthControlModel, fit!, isfitted, load_brexit

abstract type SCM end

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
    y_pre::Vector{Float64}
    y_post::Vector{Float64}
    x_pre::Matrix{Float64}
    x_post::Matrix{Float64}
    no_comps::Int64
    no_treatment_periods::Int64
    no_pretreatment_periods::Int64
    comp_labels::Vector{Union{Symbol, String}}
end


function TreatmentPanel(; data = nothing, outcome = nothing, id_var = nothing, t_var = nothing,
                        treatment::T = nothing) where T <: Pair

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

    y_pre = data[(data[!, id_var] .== treatᵢ) .&
               (data[!, t_var] .< treatₜ), outcome]
    xs_pre = data[(data[!, id_var] .!= treatᵢ) .&
               (data[!, t_var] .< treatₜ), outcome]

    y_post = data[(data[!, id_var] .== treatᵢ) .&
               (data[!, t_var] .>= treatₜ), outcome]
    xs_post = data[(data[!, id_var] .!= treatᵢ) .&
               (data[!, t_var] .>= treatₜ), outcome]

    no_treatment_periods = length(unique(data[data[!, t_var] .>= treatₜ, t_var]))
    no_pretreat_periods = length(y_pre)
    no_comps = length(unique(data[data[!, id_var] .!= treatᵢ, id_var]))

    comp_labels = unique(data[data[!, id_var] .!= treatᵢ, id_var])

    TreatmentPanel(
        data,
        outcome,
        id_var,
        t_var,
        treatment,
        y_pre,
        y_post,
        reshape(xs_pre, no_comps, no_pretreatment_periods),
        reshape(xs_post, no_comps, no_treatment_periods),
        no_comps,
        no_treatment_periods,
        no_pretreat_periods,
        comp_labels
    )
end


"""
    SynthControlModel

Takes in a `TreatmentPanel` and holds the results of calling the `fit!` method, optimal
weights for comparator units and the resulting predicted outcome of the synthetic control
unit as well as the estimated treatment effect from comparing synthetic control and observed
outcome of the treated unit.

"""

mutable struct SynthControlModel <: SCM
    treatment_panel::TreatmentPanel
    w::Vector{Float64} # Vector of weights
    ŷ::Vector{Float64} # Vector of post-treatment predicted outcomes for synthetic unit
    δ̂::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end


function SynthControlModel(tp::TreatmentPanel)
    @unpack y_pre, y_post, no_comps, no_treatment_periods = tp

    ŷ = zeros(length(y_pre) + length(y_post))
    w = zeros(no_comps)
    δ̂ = zeros(no_treatment_periods)
    p_test_res = zeros(length(ŷ), no_comps)

    SynthControlModel(tp, w, ŷ, δ̂, p_test_res)
end


"""
isfitted(s::SynthControlModel)

Check whether the `SynthControlModel` object `s` has been fitted.

# Examples
```julia-repl
julia> df = load_brexit();

julia> s = SynthControlModel(df, :country, :dateid, :realgdp, 86, "United Kingdom");

julia> isfitted(s)
false

julia> fit!(s);

julia> isfitted(s)
true

```
"""
isfitted(s::SynthControlModel) = sum(s.w) ≈ 0.0 ? false : true

"""
    fit!(s::SynthControlModel; placebo_test = false)

Fit the `SynthControlModel` `s` by finding the weights that minimize the distance between
the pre-treatment outcomes for the observational unit of interest and the weighted average
pre-treatment outcomes for unweighted units.

If `placebo_test = true` is supplied, additional placebo tests will be performed by using
every non-treated unit in the data set as the treated unit in turn and estimating the
treatment impact on this unit. Results are stored in the `p_test_res` field and can be
used as the basis for inference.
"""
function fit!(s::SynthControlModel; placebo_test = false)

    @unpack y_pre, y_post, x_pre, x_post, no_comps, no_pretreatment_periods = s.treatment_panel

    # Objective function
    # Hat tip to Mathieu Tanneau who suggested switching to a quadratic formulation
    obj(w) = dot(y_pre .- vec(w'*x_pre), y_pre .- vec(w'*x_pre)) + 1e6*(1.0 - sum(w))^2

    # Initial condition and bounds
    initial = [1/no_comps for _ in 1:no_comps]
    lower = zeros(no_comps)
    upper = ones(no_comps)

    # Get weights through optimization
    res = optimize(obj, lower, upper, initial)
    s.w = res.minimizer

    # Get estimates
    s.ŷ = vec(s.w'*[x_pre x_post])
    s.δ̂ = ([y_pre; y_post] .- s.ŷ)[no_pretreatment_periods+1:end]

    if placebo_test
        placebos = unique(s.treatment_panel.data[(s.treatment_panel.data[!,s.pid] .!= s.treat_id), s.pid])
        p = deepcopy(s)

        for n ∈ 1:p.no_comps
            # change treatment ID
            p.treat_id = placebos[n]
            # change outcome data
            p.y = p.data[p.data[!, p.pid] .== p.treat_id, p.outcome]
            p.comps = collect(reshape(p.data[(p.data[!,p.pid] .!= p.treat_id) .& (p.data[!,p.tid].<p.treat_t), p.outcome],
                    p.no_comps, p.no_pretreat_t)')
            fit!(p)
            s.p_test_res[:, n] = p.y .- p.ŷ
        end
    end
    return s
end

# Define plotting recipe
@recipe function f(s::SynthControlModel; kind = "overall")
    legend --> :topleft

    if kind == "overall"

        @series begin
            label --> s.treatment_panel.treatment[1]
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y_pre; s.treatment_panel.y_post]
        end

        @series begin
            label --> "Control"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), s.ŷ
        end

        @series begin
            label --> "Impact"
            seriestype := :bar
            seriesalpha := 0.5
            linecolor := "white"
            seriescolor := "darkgreen"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var]))[findfirst(x -> x >= s.treatment_panel.treatment[2],
                                    sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var]))):end], s.δ̂
        end

    elseif kind == "diffplot"
        @series begin
            label --> s.treatment_panel.treatment[1]
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y_pre; s.treatment_panel.y_post] .- s.ŷ
        end

        @series begin
            label --> ""
            seriestype := :scatter
            seriescolor := 1
            markerstrokecolor := "white"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y_pre; s.treatment_panel.y_post] .- s.ŷ
        end
    end
        @series begin
            label --> "Intervention"
            seriestype := :vline
            linestyle := :dash
            seriescolor := "black"
            xguide := "Time"
            yguide := "Outcome"
            [s.treatment_panel.treatment[2]]
        end

     nothing
end

# Data loading for example
"""
    load_brexit()

Returns a `DataFrame` with quarterly GDP data on 30 OECD countries, borrowed from the
analysis of the effect of Brexit on UK GDP undertaken by the [Centre for European Reform](
https://www.cer.eu/insights/cost-brexit-june-2019)
"""
function load_brexit()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"), DataFrame)
end


# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SynthControlModel)

    println(io, "\nSynthetic Control Model")
    println(io, "\tOutcome variable: ", s.treatment_panel.outcome)
    println(io, "\tTime dimension: ",string(s.treatment_panel.t_var)," with ",
            string(length(unique(s.treatment_panel.data[!,s.t_var]))), " unique values")
    println(io, "\tTreatment period: ", string(s.treatment_panel.treatment[2]))
    println(io, "\tID variable: ",string(s.treatment_panel.id_var), " with ",
            string(length(unique(s.treatment_panel.data[!,s.treatment_panel.pid]))), " unique values")
    println(io, "\tTreatment ID: ",string(s.treatment_panel.treat_id))

    if isfitted(s)
      println(io, "\tModel is fitted")
      println(io, "\tImpact estimates: ",round.(s.δ̂, digits=3))
    else
      println(io, "\n\tModel is not fitted")
    end
end

end # module
