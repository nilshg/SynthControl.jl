"""
    ADMSynthControlModel

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

mutable struct ADMSynthControlModel{T1} <: SCM
    treatment_panel::T1
    w::Vector{Float64} # Vector of weights for each non-treated observation, length N-1
    v::Matrix{Float64} # Vector of importance weights for predictors X₁, ..., Xₖ, length k 
    ŷ₁::Vector{Float64} # Vector of post-treatment predicted outcomes for synthetic unit
    τ̂::Vector{Float64} # Vector of estimated treatmet impacts
    p_test_res::Array{Float64, 2} # Matrix of outcomes for untreated under placebo
end


function ADMSynthControlModel(tp::BalancedPanel{TreatmentPanels.SingleUnitTreatment, S}) where S
    @unpack Y, N, T, W, X = tp

    ŷ₁ = zeros(T)
    w = zeros(N - 1)
    τ̂ = zeros(length_T₁(tp))
    p_test_res = zeros(N - 1, T)

    SynthControlModel{typeof(tp)}(tp, w, ŷ₁, τ̂, p_test_res)
end

function fit!(x::ADMSynthControlModel;
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

####################################################################################################
##################################     MANUAL ATTEMPT AT ADH     ###################################
####################################################################################################

using SynthControl
using Dates, DataFrames, Statistics, LinearAlgebra

df = SynthControl.load_germany()
sort!(df, [:country, :year])

transform!(groupby(df, :country),
    [:year, :gdp] => ((y, g) -> mean(skipmissing(g[1981 .<= y .<= 1990]))) => :gdp_81_90,
    [:year, :trade] => ((y, i) -> mean(skipmissing(i[1981 .<= y .<= 1990]))) => :trade_81_90,
    [:year, :infrate] => ((y, i) -> mean(skipmissing(i[1981 .<= y .<= 1990]))) => :infrate_81_90,
    [:year, :industry] => ((y, i) -> mean(skipmissing(i[1981 .<= y .<= 1990]))) => :industry_81_90, 
    [:year, :invest80] => ((y, i) -> mean(skipmissing(i[1980 .<= y .<= 1985]))) => :invest80_80_85,
    [:year, :schooling] => ((y, i) -> mean(skipmissing(i[1980 .<= y .<= 1985]))) => :schooling_80_85)

covar_list = [:gdp_81_90,  :trade_81_90, :infrate_81_90, :industry_81_90, :schooling_80_85, :invest80_80_85]

X₁_df = combine(first, groupby(df[(df.country .== "West Germany"), [covar_list; :country]], :country))
X₁ = [x for x ∈ X₁_df[1, 2:end]]

X₀_df = combine(first, groupby(df[(df.country .!= "West Germany"), [covar_list; :country]], :country, sort = true))
X₀ = Matrix(X₀_df[:, 2:end])'

w_equal = [1/size(X₀, 2) for _ ∈ 1:size(X₀, 2)]

# Weights shown in ADH ()
w_adh = [0, 0.42, 0, 0, 0, 0, 0, 0.16, 0.09, 0, 0, 0, 0, 0.11, 0, 0.22]

sqrt((X₁ .- X₀*w_adh)'*(X₁ .- X₀*w_adh))
norm(X₁ .- X₀*w_adh)
norm(X₁ .- X₀*w_equal)

V = zeros(6, 6)
V[diagind(V)] .= 1.0

sqrt((X₁ .- X₀*w_equal)'*V*(X₁ .- X₀*w_equal))
sqrt((X₁ .- X₀*w_adh)'*V*(X₁ .- X₀*w_adh))

obj1(w) = (X₁ .- X₀*w)'*(X₁ .- X₀*w)
obj2(w) = (X₁ .- X₀*w)'*(X₁ .- X₀*w) + 1e6 * (1.0 - sum(w))^2
obj3(w) = (X₁ .- X₀*w)'*(X₁ .- X₀*w) + 1e6 * (1.0 - sum(w))^2 + 1_000 * norm(w)

using Optim

initial = w_equal
upper = fill(1.0, length(w_equal))
lower = fill(0.0, length(w_equal))

w_simple = optimize(obj1, lower, upper, initial).minimizer
w_sum_constraint = optimize(obj2, lower, upper, initial).minimizer
w_sum_ridge_constraint = optimize(obj3, lower, upper, initial).minimizer

summary = DataFrame(Country = sort(unique(df[df.country .!= "West Germany", :].country)), 
            ADH_weights = w_adh,
            w_simple = round.(w_simple, digits = 2),
            w_sum_constraint = round.(w_sum_constraint, digits = 2),
            w_sum_ridge_constraint = round.(w_sum_ridge_constraint, digits = 2))

V_adh = zeros(6, 6)
covar_list = [:gdp_81_90,  :trade_81_90, :infrate_81_90, :industry_81_90, :schooling_80_85, :invest80_80_85]
V_adh[diagind(V_adh)] .= [0.442, 0.134, 0.072, 0.001, 0.107, 0.245]

sqrt((X₁.- X₀*w_adh)'*V_adh*(X₁ .- X₀*w_adh))
