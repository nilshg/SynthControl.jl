
abstract type SCM end


"""
    SimpleSCM

Takes in a `TreatmentPanel` and holds the results of calling the `fit!` method, optimal
weights for comparator units and the resulting predicted outcome of the synthetic control
unit as well as the estimated treatment effect from comparing synthetic control and observed
outcome of the treated unit.

"""
mutable struct SimpleSCM{T1} <: SCM
    treatment_panel::T1
    ω::Vector{Float64} # Vector of weights
    ŷ₁::Vector{Float64} # Vector of post-treatment predicted outcomes for synthetic unit
    τ̂::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end


function SimpleSCM(tp::BalancedPanel{SingleUnitTreatment{Continuous}})
    N = size(tp.Y, 1)
    T = size(tp.Y, 2)

    ŷ₁ = zeros(T)
    ω = zeros(N - 1)
    τ̂ = zeros(length_T₁(tp))
    p_test_res = zeros(N - 1, length_T₁(tp))

    SimpleSCM{typeof(tp)}(tp, ω, ŷ₁, τ̂, p_test_res)
end


"""
isfitted(s::SimpleSCM)

Check whether the `SimpleSCM` object `s` has been fitted.

# Examples
```julia-repl
julia> df = load_brexit();

julia> s = SimpleSCM(df, :country, :dateid, :realgdp, 86, "United Kingdom");

julia> isfitted(s)
false

julia> fit!(s);

julia> isfitted(s)
true

```
"""
isfitted(s::SimpleSCM) = sum(s.ω) ≈ 0.0 ? false : true

"""
    fit!(s::SimpleSCM; placebo_test = false)

Fit the `SimpleSCM` `s` by finding the weights that minimize the distance between
the pre-treatment outcomes for the observational unit of interest and the weighted average
pre-treatment outcomes for unweighted units.

If `placebo_test = true` is supplied, additional placebo tests will be performed by using
every non-treated unit in the data set as the treated unit in turn and estimating the
treatment impact on this unit. Results are stored in the `p_test_res` field and can be
used as the basis for inference.
"""
function fit!(s::SimpleSCM; placebo_test = false, 
    optimizer = HiGHS.Optimizer, silent = true)

    # Make sure that we have a single continuous treatment
    s isa SimpleSCM{BalancedPanel{SingleUnitTreatment{Continuous}}} ||
        throw("SimpleSCM currently only supports a single continuously treated unit")

    #!# TO DO: Check for covariates
    #isnothing(s.treatment_panel.predictors) || "Predictors currently not implemented"
    
    tp = s.treatment_panel
    
    #!# Consider how this has to change for single vs multi treatments
    T₀ = length_T₀(tp)
    y₁₀, y₁₁, yⱼ₀, yⱼ₁ = decompose_y(tp)
    J = size(yⱼ₀, 1)

    # Objective function
    # Hat tip to Mathieu Tanneau who suggested switching to a quadratic formulation and solving 
    # this with JuMP which gave a 100x speedup over Optim.jl
    # !# Check why yⱼ₀ transpose is needed - is this consistent with mathematical formulation of
    # model? Appears so, the notation for y₁₀ based on the outcome matrix and x₀ based on a vector 
    # of covariates are incompatible 
    model = Model(optimizer)
    silent && set_silent(model)
    @variable(model, 0 <= ω[1:J] <= 1)
    @constraint(model, sum(ω[i] for i ∈ 1:J) == 1)
    @objective(model, Min, (y₁₀ .- yⱼ₀'ω)'*(y₁₀ .- yⱼ₀'ω))
    
    # Get weights through optimization
    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) == JuMP.OPTIMAL
    s.ω = JuMP.value.(model[:ω])

    # Get estimates
    s.ŷ₁ = vec(s.ω'*[yⱼ₀ yⱼ₁])
    s.τ̂ = ([y₁₀; y₁₁] .- s.ŷ₁)[T₀+1:end]

    if placebo_test

        p_model = Model(optimizer)
        silent && set_silent(p_model)
        @variable(p_model, 0 <= ω[1:J-1] <= 1)
        @constraint(p_model, sum(ω[i] for i ∈ 1:J-1) == 1)

        for p ∈ 1:J
            p_y₁₀ = vec(@view yⱼ₀[p, :])
            p_y₁₁ = vec(@view yⱼ₁[p, :])
            p_yⱼ₀ = @view yⱼ₀[Not(p), :]
            p_yⱼ₁ = @view yⱼ₁[Not(p), :]

            @objective(p_model, Min, (p_y₁₀ .- p_yⱼ₀'ω)'*(p_y₁₀ .- p_yⱼ₀'ω))
            
            # Get weights through optimization
            JuMP.optimize!(p_model)
            @assert JuMP.termination_status(p_model) == JuMP.OPTIMAL
            p_ω = JuMP.value.(p_model[:ω])
            
            # Get estimates
            p_ŷ₁ = vec(p_ω'*[p_yⱼ₀ p_yⱼ₁])
            s.p_test_res[p, :] = ([p_y₁₀; p_y₁₁] .- p_ŷ₁)[T₀+1:end]
        end
    end

    return s
end

# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SimpleSCM)

    println(io, "Synthetic Control Model\n")
    println(io, "Treatment panel:")
    display(s.treatment_panel)
    
    if isfitted(s)
      println(io, "\tModel is fitted")
      println(io, "\tImpact estimates: ",round.(s.τ̂, digits=3))
      println(io, "\tATT: $(round(sum(s.τ̂)/length(s.τ̂), digits = 3))")
    else
      println(io, "\nModel is not fitted")
    end
end
