"""
    SyntheticDiD

Synthetic Differences-in-Differences (SDID) estimator as described in Arkhangelsky et al. (2021).
SDID is defined as 

    μ̂, α̂, β̂, τ̂ = argmin ∑ₙ∑ₜ(Yᵢₜ - μ - αᵢ - βₜ - Wᵢₜτ)²ω̂ᵢλ̂ₜ

where i = 1, 2, ..., N denotes the observation units in the panel, t = 1, ..., T denotes the time
periods, Y is an (N×T) matrix of outcomes, W is an (N×T) indicator matrix of treatment status, ω is
a unit weight and λ is a time period weight. 

The implementation follows the author's reference implementation in the R package `synthdid`.

NOTE: The implementation assumes that the outcomes and treatment matrices are sorted such that
treated units come last in both `Y` and `W`. It therefore checks whether `W` is sorted and swaps
rows in both `Y` and `W` to sort both matrices accordingly if required. 

"""
struct SyntheticDiD{T1}# <: SCM where T1
    tp::T1
    ω̂::Vector{Float64}
    λ̂::Vector{Float64}
    ŷ::Vector{Float64}
    τ̂::Vector{Float64}
    se_τ̂::Vector{Float64}
end

function SyntheticDiD(tp::BalancedPanel{SingleUnitTreatment{Continuous}})
    ω̂ = zeros(length(tp.is) - 1)
    λ̂ = zeros(length_T₀(tp))
    ŷ = zeros(size(tp.Y, 2))
    τ̂ = zeros(1)
    se_τ̂ = zeros(1)
    SyntheticDiD(tp, ω̂, λ̂, ŷ, τ̂, se_τ̂)
end

function isfitted(x::SyntheticDiD)
    all(!iszero, x.τ̂) || all(!iszero, x.se_τ̂)
end

function fit!(x::SyntheticDiD{BalancedPanel{SingleUnitTreatment{Continuous}}};
    se = nothing)
    (;Y, W) = x.tp
    
    T₀ = length_T₀(x.tp)
    N₁ = 1 # length(treated_ids(x)) - has to be one given type of BalancedPanel in this method
    N₀ = size(Y, 1) - N₁

    # Implementation of synthdid below assumes that treated unit comes last
    if !issorted(any(x) for x ∈ eachrow(W)) 
        treated_row = treated_ids(x.tp)
        swaprows!(Y, treated_row, size(Y, 1))
        swaprows!(W, treated_row, size(W, 1))
    end

    # Get estimate and store in SyntheticDiD object
    estimate = synthdid_estimate(Y, N₀, T₀)   
    x.τ̂ .= estimate.τ̂
    x.λ̂ .= estimate.λ̂
    x.ω̂ .= estimate.ω̂

    if !isnothing(se)
        if se != :placebo
            error("Only placebo standard error estimation is currently implemented")
        else
            V̂ = get_V̂_τ(x).V̂
            x.se_τ̂ .= √(V̂)
        end
    end

    # return SynhtDiD object
    return x
end


function get_V̂_τ(x::SyntheticDiD{BalancedPanel{SingleUnitTreatment{Continuous}}})
    (;Y, W) = x.tp

    Y_untreated = Y[Not(treated_ids(x.tp)), :]

    T₀ = length_T₀(x.tp)
    N₁ = 1
    N₀ = size(Y_untreated, 1) - N₁

    τ̂s = zeros(N₀)

    for i ∈ 1:N₀-1
        swaprows!(Y_untreated, i, size(Y_untreated, 1))
        τ̂s[i] = synthdid_estimate(Y_untreated, N₀, T₀).τ̂
    end

    swaprows!(Y_untreated, 1, size(Y_untreated, 1))
    τ̂s[end] = synthdid_estimate(Y_untreated, N₀, T₀).τ̂

    return (V̂ = var(τ̂s), τ̂s = τ̂s)
end

function synthdid_estimate(Y, N₀, T₀;
    ω_intercept = true, λ_intercept = true,
    max_iter = 10_000, max_iter_pre_sparsify = 100, sparsify_function = sparsify_weights!)

    @assert size(Y, 1) > N₀ && size(Y, 2) > T₀

    N1 = size(Y, 1) - N₀
    T1 = size(Y, 2) - T₀

    ## (1) Compute regularization parameter ζ
    # First difference of outcomes for control units before treatment time
    Δᵢₜ = Y[1:N₀, 2:T₀] .- Y[1:N₀, 1:T₀-1]
    Δ̄ = mean(Δᵢₜ) # Average first difference
    
    #!# Paper has the uncorrected standard deviation (dividing by N̨ₖₒ*(T₀-1) only) instead, 
    # this is taken from the reference R implementation
    σ̂ = √(sum((Δᵢₜ .- Δ̄).^2)/ (N₀*(T₀ - 1)-1))
    
    ζᵢ = (N1*T1)^(1/4) * σ̂
    ζₜ = 1e-6 * σ̂

    # Calculate calculated inputs
    N1 = (size(Y, 1) - N₀)
    T1 = (size(Y, 2) - T₀)

    σ̂ = std(diff(Y[1:N₀, 1:T₀], dims = 2))
    ηᵢ = (N1 * T1)^(1/4)
    ηₜ = 1e-6
    
    ζᵢ = ηᵢ * σ̂
    ζₜ = ηₜ * σ̂

    min_decrease = 1e-5 * σ̂

    # The following is the implementation for dim(X)[3] == 0, i.e. no covariates
    Yc = collapsed_form(Y, N₀, T₀)

    ## (2) Compute unit weights ω̂
    #!# Footnote 5 mentions 10e-6*σ̂ rather than the NT^(1/4) calculation above
    # Bounds: ω₀ ∈ ℝ, ω ∈ Ω = {ω ∈ ℝ₊ᴺ: ∑wᵢ = 1} - essentially 0 ≤ ωᵢ ≤ 1 with sum constraint and ω₀
    ω̂ = sc_weight_fw(Yc[:, 1:T₀]', ζᵢ;
        intercept = ω_intercept, min_decrease = min_decrease, max_iter = max_iter_pre_sparsify).w

    
    if !isnothing(sparsify_function)
        ω̂ = sc_weight_fw(Yc[:, 1:T₀]', ζᵢ; w = sparsify_function(ω̂),
            intercept = ω_intercept, min_decrease = min_decrease, max_iter = max_iter).w
    end

    ## (3) Compute time weights λ̂x
    λ̂ = sc_weight_fw(Yc[1:N₀, :], ζₜ; 
        intercept = λ_intercept, min_decrease = min_decrease, max_iter = max_iter_pre_sparsify).w

    if !isnothing(sparsify_function)
        λ̂ = sc_weight_fw(Yc[1:N₀, :], ζₜ; w = sparsify_function(λ̂),
            intercept = λ_intercept, min_decrease = min_decrease, max_iter = max_iter).w
    end

    ## (4) Estimate τ
    τ̂ = [-ω̂; fill(1/N1, N1)]' * Y * [-λ̂; fill(1/T1, T1)]

    return (; τ̂, λ̂, ω̂)
end

function collapsed_form(Y, N₀, T₀)
    N, T = size(Y)

    [[Y[1:N₀, 1:T₀] mean(Y[1:N₀, T₀+1:T], dims = 2)];
     [mean(Y[N₀+1:N, 1:T₀], dims = 1) mean(Y[N₀+1:N, T₀+1:T], dims = 2)]]
end


function sc_weight_fw(Y, ζ; 
    intercept = true, w = nothing, min_decrease = 1e-3, max_iter = 1_000)

    T₀ = size(Y, 2) - 1
    N₀ = size(Y, 1)

    w = isnothing(w) ? fill(1/T₀, T₀) : w
    
    if intercept
        Y .-= mean(Y, dims = 1)
    end

    t = 0
    vals = fill(NaN, max_iter)
    A = Y[:, 1:T₀]
    b = Y[:, T₀+1]
    η = N₀ * ζ^2

    ζ² = ζ^2

    while t < max_iter && (t < 2 || vals[t-1] - vals[t] > min_decrease^2)
        t += 1
        wₚ = fw_step(A, w, b, η)
        w = wₚ
        #!# TO DO: this just seems to be ϵ = A * w - b, figure out whether there are cases where 
        # it might not be, otherwise go for that simpler way of writing this
        #ϵ = @view(Y[1:N₀, :]) * [w; -1]
        ϵ = A * w - b
        vals[t] = ζ² * sum(w.^2) + sum(ϵ.^2) / N₀
    end

    return (w = w, vals = vals)
end


function fw_step(A, x, b, η)
    Ax = A*x # 575ns + 3 allocs, 5.3k
    half_grad = vec((Ax .- b)' * A + (η .* x)') # This is 250ns, 4 allocs
    
    i = last(findmin(half_grad)) # basically free

    d_x = -x # <50ns, but 1 allocation
    d_x[i] = 1 - x[i] # basically free
    
    all(==(0), d_x) && return x # basically free

    d_ϵ = A[:, i] .- Ax # 50ns, 2 allocations, 0.6k
    step = (-half_grad' * d_x) / (sum(d_ϵ.^2) + η * sum(d_x.^2)) # 100ns, 3 allocs, 0.8k
    constrained_step = clamp(step, 0, 1) # basically free

    return x .+ constrained_step .* d_x
end


function fw_step(A, x, b, η, α)
    Ax = A*x
    half_grad = vec((Ax .- b)' * A + (η .* x)')
    
    i = last(findmin(half_grad))

    x *= (1-α)
    x[i] += α
        
    return x
end


function swaprows!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = axes(X,2)
        X[i, k], X[j, k] = X[j, k], X[i, k]
    end
end


function sparsify_weights!(v)
    v .= ifelse.(v .<= maximum(v)/4, 0.0, v)
    v ./= sum(v)
end


# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SyntheticDiD)

    println(io, "Synthetic Difference-in-Differences Model")
    
    if isfitted(s)
      println(io, "\tModel is fitted")
      println(io, "\tImpact estimate: ",round.(only(s.τ̂), digits=3))
      !iszero(s.se_τ̂) && println(io, "\t                  (",round.(only(s.se_τ̂), digits = 2),")")
    else
      println(io, "\nModel is not fitted")
    end

    #println(io, "Treatment panel:")
    #display(s.tp)
end
