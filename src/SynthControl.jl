module SynthControl

using DataFrames, Optim, LinearAlgebra, RecipesBase
import CSV

export SynthControlModel, fit!, isfitted, load_brexit

mutable struct SynthControlModel
    data::DataFrame
    pid::Symbol
    tid::Symbol
    outcome::Symbol
    treat_t
    treat_id
    y::Vector{Float64}
    no_comps::Int
    no_pretreat_t::Int
    no_treat_t::Int
    comps::Array{Float64, 2}
    w::Vector{Float64}
    ŷ::Vector{Float64}
    δ::Vector{Float64}
    p_test_res::Array{Float64, 2}
end

SynthControlModel(data::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
    treat_t, treat_id) = SynthControlModel(
    data,
    pid,
    tid,
    outcome,
    treat_t,
    treat_id,
    data[data[pid].==treat_id, outcome], # y
    length(unique(data[pid])) - 1, # no_comps
    length(unique(data[data[tid] .< treat_t, tid])), # no_pretreat_t
    length(unique(data[data[tid] .>= treat_t, tid])), # no_treat_t
    collect(reshape(data[.!(data[pid] .== treat_id) .&( data[tid].<treat_t), outcome],
            length(unique(data[pid])) - 1, length(unique(data[data[tid] .< treat_t, tid])))'), # comps
    zeros(length(unique(data[data[tid].>=treat_t, tid]))), # ŷ (uninitialized)
    zeros(length(unique(data[pid])) - 1), # w
    zeros(length(unique(data[data[tid] .>= treat_t, tid]))), # delta
    zeros(length(unique(data[data[tid] .>= treat_t, tid])),
          length(unique(data[pid]))) # p_test_res
    )

isfitted(s::SynthControlModel) = sum(s.w) ≈ 0.0 ? false : true

function fit!(s::SynthControlModel; placebo_test = false)

    # Estimation target
    y_pre = s.y[1:s.no_pretreat_t]

    # Objective function
    obj(w) = norm(y_pre - sum((w .* s.comps')', dims=2)) + 1_000_000*abs(1-sum(w))

    # Initial condition and bounds
    initial = [1e-5 for _ in 1:s.no_comps]
    lower = zeros(s.no_comps)
    upper = ones(s.no_comps)

    # Get weights through optimization
    res = optimize(obj, lower, upper, initial)
    s.w = res.minimizer

    # Get estimates
    s.ŷ = vec(sum((s.w .* reshape(s.data[.!(s.data[s.pid] .== s.treat_id), s.outcome],
                    s.no_comps, length(unique(s.data[s.tid]))))', dims = 2))

    s.δ = (s.y .- s.ŷ)[s.no_pretreat_t+1:end]

    if placebo_test
        placebos = unique(s.data[.!(s.data[s.pid] .== s.treat_id), s.pid])
        p = deepcopy(s)

        for n ∈ 1:no_comps
            # change treatment ID
            p.treat_id = placebos[n]
            # change outcome data
            p.y = p.data[p.data[p.pid].==p.treat_id, p.outcome]
            p.comps = collect(reshape(p.data[.!(p.data[p.pid] .== p.treat_id
                                               ) .& (p.data[p.tid] .< p.treat_t), p.outcome],
                                      p.no_comps, p.no_pretreat_t)')
            fit!(p)
            p_test_res[n, :] = p.ŷ[s.no_pretreat_t+1:end]
        end
    end

    return s
end

#function plot_estimates(s::SynthControlModel)
#    plot(s.y, label = s.treat_id, legend = :topleft)
#    plot!(s.ŷ, label = "Control")
#    vline!([s.no_pretreat_t], label = "Treatment")
#    bar!(s.no_pretreat_t+1:length(s.y), s.δ, color = "cornflowerblue", alpha = 0.8,
#            label = "Estimated impact")
#end

# Define plotting recipe

@recipe function f(s::SynthControlModel)
    legend --> :topleft

    @series begin
        label --> s.treat_id
        s.y
    end

    @series begin
        label --> "Control"
        s.ŷ
    end

    @series begin
        label --> "Impact"
        seriestype := :bar
        #@show collect(1:s.no_treat_t), s.δ
        collect(s.no_pretreat_t+1:length(s.y)), s.δ
    end

     ()
end

# Data loading for example
function load_brexit()
    df = CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"))
    df = dropmissing(df[:, 1:4], disallowmissing=true)
end


end # module
