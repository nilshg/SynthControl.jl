module SynthControl

using Dates, DataFrames, Optim, LinearAlgebra, RecipesBase, Parameters
import CSV

export TreatmentPanel, SynthControlModel, fit!, isfitted, load_brexit, load_germany, load_basque, load_smoking

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
    ŷ₁::Vector{Float64} # Vector of post-treatment predicted outcomes for synthetic unit
    τ̂::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end


function SynthControlModel(tp::TreatmentPanel)
    @unpack y₁₀, y₁₁, J, T₁ = tp

    ŷ₁ = zeros(length(y₁₀) + length(y₁₁))
    w = zeros(J)
    τ̂ = zeros(T₁)
    p_test_res = zeros(length(ŷ₁), J)

    SynthControlModel(tp, w, ŷ₁, τ̂, p_test_res)
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

    isnothing(s.treatment_panel.predictors) || "Predictors currently not implemented"

    @unpack y₁₀, y₁₁, yⱼ₀, yⱼ₁, J, T₀ = s.treatment_panel

    # Objective function
    # Hat tip to Mathieu Tanneau who suggested switching to a quadratic formulation
    obj(w) = dot(y₁₀ .- vec(w'*yⱼ₀), y₁₀ .- vec(w'*yⱼ₀)) + 1e6*(1.0 - sum(w))^2

    # Initial condition and bounds
    initial = [1/J for _ in 1:J]
    lower = zeros(J)
    upper = ones(J)

    # Get weights through optimization
    res = optimize(obj, lower, upper, initial)
    s.w = res.minimizer

    # Get estimates
    s.ŷ₁ = vec(s.w'*[yⱼ₀ yⱼ₁])
    s.τ̂ = ([y₁₀; y₁₁] .- s.ŷ₁)[T₀+1:end]

    if placebo_test
        placebos = unique(s.treatment_panel.data[(s.treatment_panel.data[!,s.pid] .!= s.treat_id), s.pid])
        p = deepcopy(s)

        for n ∈ 1:p.J
            # change treatment ID
            p.treat_id = placebos[n]
            # change outcome data
            p.y = p.data[p.data[!, p.pid] .== p.treat_id, p.outcome]
            p.comps = collect(reshape(p.data[(p.data[!,p.pid] .!= p.treat_id) .& (p.data[!,p.tid].<p.treat_t), p.outcome],
                    p.J, p.no_pretreat_t)')
            fit!(p)
            s.p_test_res[:, n] = p.y .- p.ŷ₁
        end
    end
    return s
end

# Define plotting recipe
@recipe function f(s::SynthControlModel; kind = "overall")
    legend --> :topleft

    if kind == "weights"

        wp = sort(DataFrame(comp = s.treatment_panel.comp_labels,
                  weight = s.w), :weight, rev = true)

        @series begin
            seriestype := :bar
            label --> ""
            xguide --> "Donor"
            yguide --> "Weight"
            linecolor := "white"
            wp[wp.weight .> 0.05, :comp], wp[wp.weight .> 0.05, :weight]
        end

    elseif kind == "overall"

        @series begin
            label --> s.treatment_panel.treatment[1]
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y₁₀; s.treatment_panel.y₁₁]
        end

        @series begin
            label --> "Control"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), s.ŷ₁
        end

        @series begin
            label --> "Impact"
            seriestype --> :bar
            seriesalpha --> 0.5
            linecolor --> "white"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var]))[findfirst(x -> x >= s.treatment_panel.treatment[2],
                                    sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var]))):end], s.τ̂
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

    elseif kind == "diffplot"
        @series begin
            label --> s.treatment_panel.treatment[1]
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y₁₀; s.treatment_panel.y₁₁] .- s.ŷ₁
        end

        @series begin
            label --> ""
            seriestype := :scatter
            seriescolor := 1
            markerstrokecolor := "white"
            sort!(unique(s.treatment_panel.data[!, s.treatment_panel.t_var])), [s.treatment_panel.y₁₀; s.treatment_panel.y₁₁] .- s.ŷ₁
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
    end

     nothing
end

# Data loading for example
"""
    load_brexit()

Returns a `DataFrame` with quarterly GDP data on 30 OECD countries, borrowed from the
analysis of the effect of Brexit on UK GDP undertaken by the [Centre for European Reform](
https://www.cer.eu/insights/cost-brexit-june-2019)

The `TreatmentPanel` for Brexit should be specified as
```
TreatmentPanel("United Kingdom" => Date(2016, 7, 1), df;
                       outcome = :realgdp, id_var = :country, t_var = :quarter)
```
"""
function load_brexit()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"), DataFrame)
end

function load_germany()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","germany.csv"), DataFrame)
end

function load_basque()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","basque.csv"), DataFrame)
end


function load_smoking()
    CSV.read(joinpath(dirname(@__FILE__),"..","data","smoking.csv"), DataFrame)
end



# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SynthControlModel)

    println(io, "\nSynthetic Control Model")
    println(io, "\tOutcome variable: ", s.treatment_panel.outcome)
    println(io, "\tTime dimension: ",string(s.treatment_panel.t_var)," with ",
            string(length(unique(s.treatment_panel.data[!,s.treatment_panel.t_var]))), " unique values")
    println(io, "\tTreatment period: ", string(s.treatment_panel.treatment[2]))
    println(io, "\tID variable: ",string(s.treatment_panel.id_var), " with ",
            string(length(unique(s.treatment_panel.data[!,s.treatment_panel.id_var]))), " unique values")
    println(io, "\tTreatment ID: ",string(s.treatment_panel.treatment[1]))

    if isfitted(s)
      println(io, "\tModel is fitted")
      println(io, "\tImpact estimates: ",round.(s.τ̂, digits=3))
    else
      println(io, "\n\tModel is not fitted")
    end
end

end # module
