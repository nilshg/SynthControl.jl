"""
    ADHSynthControlModel

Synthetic control model as described in [Abadie, Diamond, Hainmueller
(2010)](https://economics.mit.edu/files/11859)
The implementation follows the exposition in [Abadie (2021)](https://economics.mit.edu/files/17847),
with the weight matrix V governing weights of covariates X in the fitting procedure chosen by
cross-validation based on a train/test split of the pre-treatment outcomes

Takes in a `TreatmentPanel` and holds the results of calling the `fit!` method, optimal
weights for comparator units and the resulting predicted outcome of the synthetic control
unit as well as the estimated treatment effect from comparing synthetic control and observed
outcome of the treated unit.

"""

mutable struct ADHSynthControlModel{T1} <: SCM
    treatment_panel::T1
    w::Vector{Float64} # Vector of weights for each non-treated observation, length N-1
    v::Matrix{Float64} # Vector of importance weights for predictors X₁, ..., Xₖ, length k 
    ŷ₁::Vector{Float64} # Vector of post-treatment predicted outcomes for synthetic unit
    τ̂::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end


function ADHSynthControlModel(tp::BalancedPanel{TreatmentPanels.SingleUnitTreatment, S}) where S
    @unpack Y, N, T, W, X = tp

    ŷ₁ = zeros(T)
    w = zeros(N - 1)
    τ̂ = zeros(length_T₁(tp))
    p_test_res = zeros(N - 1, T)

    SynthControlModel{typeof(tp)}(tp, w, ŷ₁, τ̂, p_test_res)
end

function fit!(x::ADHSynthControlModel;
    train_share = 0.7)

    k = length(X)

    # Length of train/test period in pre-intervention data
    T₀_train = round(Int, train_share * length_T₀(x))
    T₀_test = round(Int, (1 - train_share) * length_T₀(x))

    # Weight matrix v is size (k × k), simplest choice for startes is inverse variance of X
    v = zeros(k, k)
    v[diagind(v)] .= [1/var(xᵢ) for x ∈ eachcol(x)]
    
    # Calculate weights on training data
    #!# Closure that captures v?
    obj(w) = (x₁ .- x₀*w)' * v * (x₁ .- x₀*w) + 1e6*(1.0 - sum(w))^2

    # Target is minimizing the distance between predicted and actual outcome for treated unit 
    # in test period
    function target(v_diag)

        v[diagind(v)] .= v_diag

        # Get weights given proposal for v
        w′ = optimize(obj, initial, upper, lower).minimizer

        # Calculate loss in test period
        (y₁_test .- y₀_test * w′)*(y₁_test .- y₀_test * w′)
    end 

    v = optimize(target, initial, lower, upper)
end