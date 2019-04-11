# To consider
# * Multiple treatment units
# * Exclude certain unit from comparator set
# * StatsModels API
module SynthControl

using DataFrames, Optim, LinearAlgebra, Plots

export fit!, plot_estimates, SynthControlModel

mutable struct SynthControlModel
    data::DataFrame
    pid::Symbol # group identifier
    tid::Symbol # time identifier
    outcome::Symbol # variable of interest identifier
    treat_id # identifier for treated observations
    treat_t # first period of treatment
    no_treatment_t::Int64
    no_pretreatment_t::Int64
    no_comps::Int64
    y::Vector{Float64} # outcome observations
    ŷ::Vector{Float64} # predictions
    w::Vector{Float64} # weights
    p_test::Array{Float64, 2} # placebo tests
end

isfitted(s::SynthControlModel) = sum(s.w) ≈ 0.0 ? false : true

SynthControlModel(df::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
                    treat_id, treat_t) = SynthControlModel(df::DataFrame,
                    pid::Symbol,
                    tid::Symbol,
                    outcome::Symbol,
                    treat_id,
                    treat_t,
                    0,
                    0,
                    0,
                    df[df[pid] .== treat_id, outcome], # y - vector of outcomes
                    zeros(length(unique(df[df[tid].>=treat_t,tid]))), # ŷ - predictions
                    zeros(length(unique(df[.!(df[pid] .== treat_id), pid]))), # w - weights
                    zeros(length(unique(df[df[tid] .>= treat_t, tid])),
                          length(unique(df[.!(df[pid] .== treat_id), pid])))
                    )

function get_preds_weights(s::SynthControlModel; placebo = "")

    treat_id = isempty(placebo) ? s.treat_id : placebo

    df = s.data

    # Get dimensions of problem
    s.no_pretreatment_t = length(unique(df[df[s.tid] .< s.treat_t, s.tid]))
    s.no_treatment_t, s.no_comps = size(s.p_test)

    # Create comparator matrix
    comps = df[.!(df[s.pid] .== treat_id) .& (df[s.tid] .< s.treat_t), s.outcome]
    comps = collect(reshape(comps, s.no_comps, s.no_pretreatment_t)')

    # Pre treatment outcomes
    y₀ = s.y[1:s.no_pretreatment_t]

    # Objective function
    obj(w) = norm(y₀ - sum((w .* comps')', dims = 2)) + 1_000_000*abs(1-sum(w))

    # Initial condition and bounds
    initial = [0.9 for _ in 1:s.no_comps]
    lower = zeros(s.no_comps)
    upper = ones(s.no_comps)

    # Get weights through optimization
    res = optimize(obj, lower, upper, zeros(s.no_comps))
    w = res.minimizer

    # Get predictions
    ŷ = vec(sum((w .* reshape(df[.!(df[s.pid] .== treat_id), s.outcome],
                    s.no_comps, length(unique(df[s.tid]))))', dims = 2))

    return w, ŷ
end

function fit!(s::SynthControlModel; placebo_test = true)

    s.w, s.ŷ = get_preds_weights(s)

    if placebo_test
        # Standard Errors
        placebos = unique(df[.!(df[s.pid] .== s.treat_id), s.pid])
        for (i, p) ∈ enumerate(placebos)
            ŷ = get_preds_weights(s, placebo = p)[2]
            s.p_test[:, i] = ŷ[s.no_pretreatment_t+1:end]
        end
    end

    s
end

function plot_estimates(s::SynthControlModel)
    plot(s.y, label = s.treat_id, legend = :topleft)
    plot!(s.ŷ, label = "Control")
    vline!([s.no_pretreatment_t], label = "Treatment period")
end


end # module
