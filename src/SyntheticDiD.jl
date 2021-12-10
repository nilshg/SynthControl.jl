"""
    SyntheticDiD

Synthetic Differences-in-Differences (SDID) estimator as described in Arkhangelsky et al. (2021).
SDID is defined as 

    μ̂, α̂, β̂, τ̂ = argmin ∑ₙ∑ₜ(Yᵢₜ - μ - αᵢ - βₜ - Wᵢₜτ)²ω̂ᵢλ̂ₜ

where i = 1, 2, ..., N denotes the observation units in the panel, t = 1, ..., T denotes the time
periods, Y is an (N×T) matrix of outcomes, W is an (N×T) indicator matrix of treatment status, ω is
a unit weight and λ is a time period weight. 

"""
struct SyntheticDiD{T1}# <: SCM where T1
    tp::T1
    ω̂::Vector{Float64}
    λ̂::Vector{Float64}
    ŷ::Vector{Float64}
    τ̂::Vector{Float64}
    se_τ̂::Vector{Float64}
end

function SyntheticDiD(x::BalancedPanel{SingleUnitTreatment{Continuous}})
    ω̂ = zeros(length(x.is) - 1)
    λ̂ = zeros(length_T₀(x))
    ŷ = zeros(size(x.Y, 2))
    τ̂ = zeros(1)
    se_τ̂ = zeros(1)
    SyntheticDiD(x, ω̂, λ̂, ŷ, τ̂, se_τ̂)
end

function fit!(x::SyntheticDiD{BalancedPanel{SingleUnitTreatment{Continuous}}})

    (;Y, W) = x.tp
    y₁₀, y₁₁, y₀₀, y₀₁ = decompose_y(x.tp)
    T₀ = length_T₀(x.tp)
    T₁ = length_T₁(x.tp)
    Nₜᵣ = 1 # length(treated_ids(x)) - has to be one given type of BalancedPanel in this method
    Nₖₒ = size(Y, 1) - Nₜᵣ

    ### Algorithm presented in Arkhangelsky et al. (2021), v4 arXiv

    ## (1) Compute regularization parameter ζ
    # First difference of outcomes for control units before treatment time
    Δᵢₜ = Y[Not(treated_ids(x.tp)), 2:T₀] .- Y[Not(treated_ids(x.tp)), 1:T₀-1]
    Δ̄ = mean(Δᵢₜ) # Average first difference

    σ̂ = √(sum((Δᵢₜ .- Δ̄)'*(Δᵢₜ .- Δ̄)) / (Nₖₒ*(T₀ - 1)))

    ζ = (Nₜᵣ*T₁)^(1/4) * σ̂


    ## (2) Compute unit weights ω̂
    #!# Note - paper seems to generalise to multiple treated units already, here using one unit for
    # now but should be 1/Nₖₒ * sum(Y[treated_ids(x), t]) at each time period
    #!# Weights produced here are very non-sparse due to ζ*norm regularization - double check!
    function ℓᵢ(ω) 
        ω₀ = ω[1]; ω = @view ω[2:end]
        (ω₀ .+ y₀₀'*ω .- y₁₀)'*(ω₀ .+ y₀₀'*ω .- y₁₀) + 1e12*(1.0 .- sum(ω))^2 + ζ^2*T₀*norm(ω)
    end

    # Bounds: ω₀ ∈ ℝ, ω ∈ Ω = {ω ∈ ℝ₊ᴺ: ∑wᵢ = 1} - essentially 0 ≤ ωᵢ ≤ 1 with sum constraint and ω₀
    # unbounded
    #!# again this will differ for multiple treated units 
    lower = [-Inf; fill(0.0, Nₖₒ)]
    upper = [Inf; fill(1.0, Nₖₒ)]
    initial = [0.0; [1/Nₖₒ for _ ∈ 1:Nₖₒ]]

    res = optimize(ℓᵢ, lower, upper, initial)

    ω̂₀ = res.minimizer[1]  #!# this is probably actually unused?
    ω̂ = res.minimizer[2:end]
    x.ω̂ .= ω̂
    
    ## (3) Compute time weights λ̂
    # Note small regularization added as per footnote 3 of the paper
    function ℓₜ(λ)
        λ₀ = λ[1]; λ = @view λ[2:end]
        ȳ₀₁ = vec(sum(y₀₁, dims = 2) ./ size(y₀₁, 2)) # Average post-treatment outcome for control
        (λ₀ .+ y₀₀*λ .- ȳ₀₁)'*(λ₀ .+ y₀₀*λ .- ȳ₀₁) + 1e12*(1.0 .- sum(λ))^2 + ζ^2*Nₖₒ*norm(λ)
    end
    
    lower = [-Inf; fill(0.0, T₀)]
    upper = [Inf; fill(1.0, T₀)]
    initial = [0.0; [1/T₀ for _ ∈ 1:T₀]]

    res = optimize(ℓₜ, lower, upper, initial)

    λ̂₀ = res.minimizer[1]  #!# this is probably actually unused?
    λ̂ = res.minimizer[2:end]

    x.λ̂ .= λ̂
    
    ## (4) Run regression
    # Need to estimate fixed effect model with time/unit fixed effects and 
    #!# To do: for covariates we want to partial them out before this regression
    # Add unit weights - 1 for treated unit, ω̂ᵢ for untreated
    ω_df = DataFrame(Symbol(x.tp.id_var) => [treated_labels(x.tp); x.tp.is[x.tp.is .!= treated_labels(x.tp)]], :ω̂ => [1.0; ω̂])
    df = leftjoin(x.tp.df, ω_df, on = x.tp.id_var)   

    # Add time weights 
    λ_df = DataFrame(Symbol(x.tp.t_var) => x.tp.ts, :λ̂ => [fill(1/T₁, T₁); λ̂])
    df = leftjoin(df, λ_df, on = x.tp.t_var)

    # Total weights - ωλ
    df.sdid_weight = df.ω̂ .* df.λ̂

    # Add treatment dummy based on treatment matrix
    df.W = reshape(x.tp.W', length(x.tp.ts)*length(x.tp.is))

    sdid_reg = reg(df,
        term(x.tp.outcome_var) ~ term("W") + fe(term(x.tp.id_var)) + fe(term(x.tp.t_var)),
        weights = :sdid_weight, Vcov.robust())

    # Treatment effect is estimated as the coefficient on treatment dummy
    τ̂ = coef(sdid_reg)
    x.τ̂ .= τ̂
    x.se_τ̂ .= stderror(sdid_reg)

    # Counterfactual outcome is actual outcome plus treatment effect
    ŷ₁₁ = y₁₁ .+ τ̂
    x.ŷ .= [y₁₀; ŷ₁₁]
end


#!# TO DO - show method for this type