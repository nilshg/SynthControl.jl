
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

    ŷ₁ = zeros(T)
    w = zeros(N - 1)
    τ̂ = zeros(length_T₁(tp))
    p_test_res = zeros(N - 1, T)

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

    # Make sure that we have a single continuous treatment
    s isa SynthControlModel{BalancedPanel{TreatmentPanels.SingleUnitTreatment, TreatmentPanels.ContinuousTreatment}} ||
        throw("fit! currently only supports a single continuously treated unit")

    #!# TO DO: Check for covariates
    #isnothing(s.treatment_panel.predictors) || "Predictors currently not implemented"
    
    tp = s.treatment_panel
    @unpack N, T, ts, is = tp
    #!# Consider how this has to change for single vs multi treatments
    T₀ = length_T₀(tp)
    y₁₀, y₁₁, yⱼ₀, yⱼ₁ = decompose_y(tp)
    J = tp.N - 1

    # Objective function
    # Hat tip to Mathieu Tanneau who suggested switching to a quadratic formulation
    # !# Check why yⱼ₀ transpose is needed - is this consistent with mathematical formulation of
    # model? Appears so, the notation for y₁₀ based on the outcome matrix and x₀ based on a vector 
    # of covariates are incompatible 
    obj(w) = (y₁₀ .- yⱼ₀'w)'*(y₁₀ .- yⱼ₀'w) + 1e6(1.0 -sum(w))^2

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

        initial = [1/(J-1) for _ in 1:(J-1)]
        lower = zeros(J-1)
        upper = ones(J-1)

        for p ∈ 1:J
            p_y₁₀ = vec(@view yⱼ₀[p, :])
            p_y₁₁ = vec(@view yⱼ₁[p, :])
            p_yⱼ₀ = @view yⱼ₀[Not(p), :]
            p_yⱼ₁ = @view yⱼ₁[Not(p), :]
            @show size.([p_y₁₀, p_y₁₁, p_yⱼ₀, p_yⱼ₁])

            p_obj(w) = dot(p_y₁₀ .- vec(w'*p_yⱼ₀), p_y₁₀ .- vec(w'*p_yⱼ₀)) + 1e6*(1.0 - sum(w))^2
            
            # Get weights through optimization
            res = optimize(p_obj, lower, upper, initial)
            p_w = res.minimizer
            @show size(p_w), size([p_yⱼ₀ p_yⱼ₁])
            
            # Get estimates
            p_ŷ₁ = vec(p_w'*[p_yⱼ₀ p_yⱼ₁])
            s.p_test_res[p, :] = [p_y₁₀; p_y₁₁] .- p_ŷ₁
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
