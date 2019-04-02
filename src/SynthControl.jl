# To consider
# * Multilple treatment units
# * Exclude certain unit from comparator set
# * StatsModels API

#To delete - test data from CER Brexit study
df = CSV.read("C:/Users/LONNG9/Desktop/Input data - cost of Brexit 2018q3.csv")

module SynthControl

using CSV, DataFrames, Optim, LinearAlgebra, Plots


mutable struct SynthControlModel
    data::DataFrame
    pid::Symbol
    tid::Symbol
    outcome::Symbol
    treat_id
    treat_t
    y::Vector{Float64} # outcome observations
    ŷ::Vector{Float64} #
    w::Vector{Float64}
    p_test_rest::Array{Float64, 2}
end



isfitted(m::SynthControlModel) = sum(m.w) ≈ 0.0 ? false : true

SynthControlModel(df::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
                    treat_id, treat_t) = SynthControlModel(df::DataFrame,
    pid::Symbol, tid::Symbol, outcome::Symbol, treat_id, treat_t,
    df[df[pid] .== treat_id, outcome], zeros(length(unique(df[df[tid].>=treat_t,tid]))),
    zeros(length(unique(df[.!(df[pid].==treat_id),pid]))))


s = SynthControlModel(df, :country, :dateid, :realgdp, "United Kingdom", 86)

function fit!(s::SynthControlModel)
    df = s.data

    # Get dimensions of problem
    no_comps = length(unique(df[s.pid])) - 1
    no_pretreat_t = length(unique(df[df[s.tid] .< s.treat_t, s.tid]))
    no_treat_t = length(unique(df[df[s.tid] .>= s.treat_t, s.tid]))

    # Estimation target
    y_pre = s.y[1:no_pretreat_t]

    # Create comparator matrix
    comps = df[.!(df[s.pid] .== s.treat_id) .& (df[s.tid] .< s.treat_t), s.outcome]
    comps = collect(reshape(comps, no_comps, no_pretreat_t)')

    # Objective function
    obj(w) = norm(y_pre - sum((w .* comps')', dims=2)) + 1_000_000*abs(1-sum(w))

    # Initial condition and bounds
    initial = [1. / no_comps for _ in 1:no_comps]
    lower = zeros(no_comps)
    upper = ones(no_comps)

    # Get weights through optimization
    res = optimize(obj, lower, upper, initial)
    s.w = res.minimizer

    s.ŷ = vec(sum(s.w .* reshape(df[.!(df[s.pid] .== s.treat_id), s.outcome],
                    no_comps, length(unique(df[s.tid]))), dims = 1))
end

function estimate_model(s::SynthControlModel)
    df::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
    treat_id, treat_t)

    fit!(s)

    # Standard Errors
    placebos = unique(df[.!(df[pid] .== treat_id), pid])
    p_test_res = Array{Float64}(undef, length(placebos), no_treat_t)
    for (i, placebo) ∈ enumerate(placebos)
        p_test_res[i, :] = fit(df, pid, tid, outcome, placebo, t_treat)[1]
    end

    return δ, ŷ, w, p_test_res
end

δ, ŷ, w, placebos = estimate_model(df, pid, tid, outcome, "United Kingdom", 86)

function plot_estimates()
    plot(y, label = treat_id, legend = :topleft)
    plot!(ŷ, label = "Control")
    vline!([no_pretreat_t], label = "Treatment")
end

end # module

isfitted(s)

fit!(s)

isfitted(s)
