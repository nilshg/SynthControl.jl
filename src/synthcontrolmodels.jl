
abstract type SCM end


"""
    SynthControlModel

Takes in a `TreatmentPanel` and holds the results of calling the `fit!` method, optimal
weights for comparator units and the resulting predicted outcome of the synthetic control
unit as well as the estimated treatment effect from comparing synthetic control and observed
outcome of the treated unit.

"""

mutable struct SynthControlModel{T1} <: SCM
    treatment_panel::T1
    w::Vector{Float64} # Vector of weights
    ŷ₁::Vector{Float64} # Vector of post-treatment predicted outcomes for synthetic unit
    τ̂::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end


function SynthControlModel(tp::BalancedPanel{TreatmentPanels.SingleUnitTreatment, S}) where S
    @unpack Y, N, T, W = tp

    ŷ₁ = zeros(N + T)
    w = zeros(N - 1)
    τ̂ = zeros(TreatmentPanels.length_T₁(tp))
    p_test_res = zeros(length(ŷ₁), N - 1)

    SynthControlModel{typeof(tp)}(tp, w, ŷ₁, τ̂, p_test_res)
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

    #!# TO DO: Check for covariates
    #isnothing(s.treatment_panel.predictors) || "Predictors currently not implemented"

    @unpack Y, W, N, T, ts, is = s.treatment_panel
    T₀ = TreatmentPanels.length_T₀(s.treatment_panel)
    y₁₀ = @view s.treatment_panel.Y[TreatmentPanels.treated_ids(s.treatment_panel), 1:T₀]
    y₁₁ = @view s.treatment_panel.Y[TreatmentPanels.treated_ids(s.treatment_panel), T₀+1:end]
    yⱼ₀ = @view s.treatment_panel.Y[Not(TreatmentPanels.treated_ids(s.treatment_panel)), 1:T₀]
    yⱼ₁ = @view s.treatment_panel.Y[Not(TreatmentPanels.treated_ids(s.treatment_panel)), T₀+1:end]
    J = s.treatment_panel.N - 1

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


# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SynthControlModel)

    println(io, "\nSynthetic Control Model\n")
    println(io, "Treatment panel:")
    display(s.treatment_panel)

    if isfitted(s)
      println(io, "\tModel is fitted")
      println(io, "\tImpact estimates: ",round.(s.τ̂, digits=3))
    else
      println(io, "\nModel is not fitted")
    end
end
