module SynthControl

using DataFrames, Optim, LinearAlgebra, RecipesBase
import CSV

export SynthControlModel, fit!, isfitted, load_brexit

mutable struct SynthControlModel
    # complete dataframe
    data::DataFrame
    # column with region names
    pid::Symbol
    # column with integer timesteps
    tid::Symbol
    # column with response variable
    outcome::Symbol
    # treatment timestep (a value in tid)
    # the treatment timestep will be considered "post-treatment"
    treat_t
    # treated region (a value in pid)
    treat_id
    # response vector for the treated region
    y::Vector{Float64}
    # number of untreated regions
    no_comps::Int
    # number of pre-treatment timesteps (excludes treatment timestep)
    no_pretreat_t::Int
    # number of post-treatment timesteps (includes treatment timestep)
    no_treat_t::Int
    # matrix of pre-treatment response values for untreated regions,
    # where each row is a timestep and each column is a region
    comps::Array{Float64, 2}
    # a vector of weights for each untreated region
    w::Vector{Float64}
    # a vector of predictions for y
    ŷ::Vector{Float64}
    # a vector of impact estimates for each post-treatment timestep
    δ::Vector{Float64}
    # for placebo test
    p_test_res::Array{Float64, 2}
end

function SynthControlModel(data::DataFrame, pid::Symbol, tid::Symbol, outcome::Symbol,
    treat_t, treat_id)

    y = data[data[!, pid] .== treat_id, outcome]
    no_comps = length(unique(data[!, pid])) - 1
    no_pretreat_t = length(unique(data[data[!, tid] .< treat_t, tid]))
    no_treat_t = length(unique(data[data[!, tid] .>= treat_t, tid]))
    # vector of pre-treatment response values for untreated regions
    comp_outcomes = data[(data[!, pid] .!= treat_id) .& (data[!, tid] .< treat_t), outcome]
    comps = collect(reshape(comp_outcomes, no_comps, no_pretreat_t)')
    w = zeros(no_comps)
    ŷ = zeros(length(y))
    δ = zeros(no_treat_t)
    p_test_res = zeros(length(y), no_comps)

    SynthControlModel(data, pid, tid, outcome, treat_t, treat_id, y, no_comps,
        no_pretreat_t, no_treat_t, comps, w, ŷ, δ, p_test_res)
end

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
    s.ŷ = vec(sum((s.w .* reshape(s.data[(s.data[!,s.pid] .!= s.treat_id), s.outcome],
                    s.no_comps, length(unique(s.data[!,s.tid]))))', dims = 2))

    s.δ = (s.y .- s.ŷ)[s.no_pretreat_t+1:end]

    if placebo_test
        placebos = unique(s.data[(s.data[!,s.pid] .!= s.treat_id), s.pid])
        p = deepcopy(s)

        for n ∈ 1:p.no_comps
            # change treatment ID
            p.treat_id = placebos[n]
            # change outcome data
            p.y = p.data[p.data[!, p.pid] .== p.treat_id, p.outcome]
            p.comps = collect(reshape(p.data[(p.data[!,p.pid] .!= p.treat_id) .& (p.data[!,p.tid].<p.treat_t), p.outcome],
                    p.no_comps, p.no_pretreat_t)')
            fit!(p)
            s.p_test_res[:, n] = p.y .- p.ŷ
        end
    end
    return s
end

# Define plotting recipe
@recipe function f(s::SynthControlModel)
    legend --> :topleft

    @series begin
        label --> s.treat_id
        sort!(unique(s.data[!, s.tid])), s.y
    end

    @series begin
        label --> "Control"
        sort!(unique(s.data[!, s.tid])), s.ŷ
    end

    @series begin
        label --> "Intervention"
        seriestype := :vline
        linestyle := :dash
        [s.treat_t]
    end

    @series begin
        label --> "Impact"
        seriestype := :bar
        alpha := 0.3
        colour := "salmon"
        minimum(s.data[!, s.tid]) .+ collect(s.no_pretreat_t+1:length(s.y)), s.δ
    end

     nothing
end

# Data loading for example
function load_brexit()
    df = CSV.read(joinpath(dirname(@__FILE__),"..","data","brexit.csv"))
    df = dropmissing(df[:, 1:4], disallowmissing=true)
end

# Pretty printing
import Base.show

function show(io::IO, ::MIME"text/plain", s::SynthControlModel)

    println("Synthetic Control Model\n")
    println("Outcome variable: ",string(s.outcome))
    println("Time dimension: ",string(s.tid)," with ",
            string(length(unique(s.data[!,s.tid]))), " unique values")
    println("Treatment period: ",string(s.treat_t))
    println("ID variable: ",string(s.pid), " with ",
            string(length(unique(s.data[!,s.pid]))), " unique values")
    println("Treatment ID: ",string(s.treat_id))

    if isfitted(s)
      println("Model is fitted")
      println("Impact estimates: ",round.(s.δ, digits=3))
    else
      println("\nModel is not fitted")
    end
end

end # module
