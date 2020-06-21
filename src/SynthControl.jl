module SynthControl

using DataFrames, Optim, LinearAlgebra, RecipesBase
import CSV

export SynthControlModel, fit!, isfitted, load_brexit

"""
    SynthControlModel

Holds the data, model information, and results of a synthetic control model

# Members

- `data`: The dataset on which the model is esitmated
- `outcome`: Column name of the outcome variable of interest in the data
- `pid`: Column name of the panel id variable in the data
- `tid`: Column name of the time period variable in the data
- `treat_id`:  The name of the treated unit - should be an element of `data.pid`
- `treat_t`: The first period of treatment - should be an element of `data.tid`

Currently, only a single treatment unit and a single treatment period are supported.


"""
mutable struct SynthControlModel
    data::DataFrame # Data set
    outcome::Symbol # Variable of interest
    pid::Symbol # Column name of column holding panel id variable
    tid::Symbol # Column name of column holdign time id variable
    treat_id # Identifier of treated unit
    treat_t # First treatment period
    y::Vector{Float64} # Vector of outcomes
    no_comps::Int # Number of untreated units
    no_pretreat_t::Int # Number of pre-treatment periods
    no_treat_t::Int # Number of treatment periods
    comps::Array{Float64, 2} # Matrix of pre-treatment outcomes for untreated units
                             # dimension (no_comps × no_pretreat_t)
    w::Vector{Float64} # Vector of weights
    ŷ::Vector{Float64} # Vector of post-treatment outcomes for synthetic unit
    δ::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end

function SynthControlModel(data::DataFrame, outcome::Symbol, pid::Symbol, tid::Symbol,
    treat_id, treat_t)

    y = data[data[!, pid] .== treat_id, outcome]
    no_comps = length(unique(data[! ,pid])) - 1
    no_pretreat_t = length(unique(data[data[! ,tid] .< treat_t, tid]))
    no_treat_t = length(unique(data[data[! ,tid] .>= treat_t, tid]))
    comp_outcomes = data[(data[!,pid] .!= treat_id) .& (data[!,tid].<treat_t), outcome]
    comps = collect(reshape(comp_outcomes, no_comps, no_pretreat_t)')
    ŷ = zeros(length(y))
    w = zeros(no_comps)
    δ = zeros(no_treat_t)
    p_test_res = zeros(length(y), no_comps)

    SynthControlModel(data, outcome, pid, tid, treat_id, treat_t, y, no_comps,
        no_pretreat_t, no_treat_t, comps, w, ŷ, δ, p_test_res)
end

function SynthControlModel(data::DataFrame; outcome::Symbol = nothing,
    pid::Symbol = nothing, tid::Symbol = nothing, treat_id = nothing,
    treat_t = nothing)

    if isnothing(outcome)
        error(ArgumentError("""
        Please specify outcome, the nameof the column in your dataset holding the
        outcome variable of interest in your dataset
        """))
    elseif isnothing(pid)
        ArgumentError("""
            Please specify pid, the name of the column in your data set holding the
            identifier of your units of observation
            """)
    elseif isnothing(tid)
        error(ArgumentError("""
            Please specify tid, the name of the column in your dataset holding the
            time dimension
            """))
    elseif isnothing(treat_id)
        ArgumentError("""
        Please specify treat_id, the identifier of the treated unit(s) in your
        dataset
        """)
    elseif isnothing(treat_t)
        ArgumentError("""
            Please specify treat_t, the first treatment period
            """)
    else
        SynthControlModel(data, outcome, pid, tid, treat_id, treat_t)
    end
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

    # Estimation target
    y_pre = s.y[1:s.no_pretreat_t]

    # Objective function
    obj(w) = norm(y_pre .- s.comps*w) + 1_000_000*abs(1-sum(w))

    # Initial condition and bounds
    initial = [1e-5 for _ in 1:s.no_comps]
    lower = zeros(s.no_comps)
    upper = ones(s.no_comps)

    # Get weights through optimization
    res = optimize(obj, lower, upper, initial)
    s.w = res.minimizer

    # Get estimates
    s.ŷ = vec(sum((s.w .* reshape(s.data[(s.data[!,s.pid] .!= s.treat_id), s.outcome],
                    s.no_comps, length(unique(s.data[!,s.tid]))))', dims = 2))

    s.δ = (s.y .- s.ŷ)[s.no_pretreat_t+1:end]

    if placebo_test
        placebos = unique(s.data[(s.data[!,s.pid] .!= s.treat_id), s.pid])
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
@recipe function f(s::SynthControlModel)
    legend --> :topleft

    @series begin
        label --> s.treat_id
        sort!(unique(s.data[!, s.tid])), s.y
    end

    @series begin
        label --> "Control"
        sort!(unique(s.data[!, s.tid])), s.ŷ
    end

    @series begin
        label --> "Intervention"
        seriestype := :vline
        linestyle := :dash
        [s.treat_t]
    end

    @series begin
        label --> "Impact"
        seriestype := :bar
        alpha := 0.3
        colour := "salmon"
        minimum(s.data[!, s.tid]) .+ collect(s.no_pretreat_t+1:length(s.y)), s.δ
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
    df = DataFrame(CSV.File(joinpath(dirname(@__FILE__),"..","data","brexit.csv")))
    df = dropmissing(df[:, 1:4], disallowmissing=true)
end

# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SynthControlModel)

    println("Synthetic Control Model\n")
    println("Outcome variable: ",string(s.outcome))
    println("Time dimension: ",string(s.tid)," with ",
            string(length(unique(s.data[!,s.tid]))), " unique values")
    println("Treatment period: ",string(s.treat_t))
    println("ID variable: ",string(s.pid), " with ",
            string(length(unique(s.data[!,s.pid]))), " unique values")
    println("Treatment ID: ",string(s.treat_id))

    if isfitted(s)
      println("Model is fitted")
      println("Impact estimates: ",round.(s.δ, digits=3))
    else
      println("\nModel is not fitted")
    end
end

end # module
