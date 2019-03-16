#To delete - test data from CER Brexit study
df = CSV.read("C:/Users/Nils-Holger/Desktop/Input data - cost of Brexit 2018q3.csv")

treat_id = "United Kingdom"
pid = :country
tid = :dateid
outcome = :realgdp
t_treat = 86

module SynthControl

import Base:show
using DataFrames, Optim, LinearAlgebra, Plots
export estimate_model

function fit(df::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
    treat_id, treat_t)

    # Get dimensions of problem
    no_comps = length(unique(df[pid])) - 1
    no_pretreat_t = length(unique(df[df[tid] .< t_treat, tid]))
    no_treat_t = length(unique(df[df[tid] .>= t_treat, tid]))

    # Estimation target
    y = df[df[pid] .== treat_id, outcome]
    y_pre = y[1:no_pretreat_t]

    # Create comparator matrix
    comps = df[.!(df[pid] .== treat_id) .& (df[tid] .< t_treat), outcome]
    comps = collect(reshape(comps, no_comps, no_pretreat_t)')

    # Objective function
    obj(w) = norm(y_pre - sum((w .* comps')', dims=2)) + 1_000_000*abs(1-sum(w))

    # Initial condition and bounds
    initial = [1. / no_comps for _ in 1:no_comps]
    lower = zeros(no_comps)
    upper = ones(no_comps)

    # Get weights through optimization
    res = optimize(obj, lower, upper, initial)
    w_res = res.minimizer

    ŷ = sum((w_res .* reshape(df[.!(df.country .== "United Kingdom"), :realgdp],
                    no_comps, length(unique(df.dateid))))', dims = 2)

    return (y-ŷ)[no_pretreat_t+1:end], ŷ, w_res
end

function estimate_model(df::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
    treat_id, treat_t)

    δ, ŷ, w = fit(df, pid, tid, outcome, treat_id, treat_t)

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
