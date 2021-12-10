[![CI](https://github.com/nilshg/SynthControl.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/nilshg/SynthControl.jl/actions/workflows/ci.yml)

# SynthControl.jl
Julia package for synthetic control methods

The package is currently at alpha stage - it currently offers a basic `SyntheticControlModel`
implementation which uses all pre-treatment periods and no other covariates in constructing the
control unit. 

## Installation

The package is registered in the general registry, installation therefore works through the Pkg REPL

```
pkg> add SynthControl
```

## Usage

### Fitting a `SimpleSCM` model

The package includes example data borrowed from the [CER Brexit study](
  https://www.cer.eu/insights/cost-brexit-june-2018):

```
julia> using SynthControl, Dates

julia> df = load_brexit()
897×3 DataFrame
│ Row │ country        │ quarter    │ realgdp │
│     │ String         │ Date       │ Float64 │
├─────┼────────────────┼────────────┼─────────┤
│ 1   │ Australia      │ 2009-01-01 │ 1.04    │
│ 2   │ Austria        │ 2009-01-01 │ -1.53   │
│ 3   │ Belgium        │ 2009-01-01 │ -1.15   │
⋮
│ 894 │ Sweden         │ 2018-07-01 │ 22.48   │
│ 895 │ Switzerland    │ 2018-07-01 │ 14.35   │
│ 896 │ United Kingdom │ 2018-07-01 │ 15.72   │
│ 897 │ United States  │ 2018-07-01 │ 19.32   │
```

The package defines a `SimpleSCM` type, instances of which can be constructed from a
`TreatmentPanel` object from the package
[`TreatmentPanels`](https://github.com/nilshg/TreatmentPanels.jl). The `TreatmentPanel` is
constructed from a `DataFrame` and a specification of treatment assignment. The `SimpleSCM` model
constructs a synthetic control unit based only 

The example data set includes quarterly GDP for a number of OECD countries, and
we are interested in estimating the impact of the Brexit vote in Q2 2016 on GDP
in the UK:

```
julia> bp = BalancedPanel(df, "United Kingdom" => Date(2016, 7, 1); id_var = :country, t_var = :quarter, outcome_var = :realgdp)
Balanced Panel - single treated unit, continuous treatment
    Treated unit: United Kingdom
    Number of untreated units: 22
    First treatment period: 2016-07-01
    Number of pretreatment periods: 30
    Number of treatment periods: 9


julia> s = SimpleSCM(bp)

Synthetic Control Model

Treatment panel:
Balanced Panel - single treated unit, continuous treatment
    Treated unit: United Kingdom
    Number of untreated units: 22
    First treatment period: 2016-07-01
    Number of pretreatment periods: 30
    Number of treatment periods: 9

Model is not fitted
```

The output indicates that the model is not fitted, that is we have at this stage
only defined the basic model structure. We can fit the model using the `fit!`
function, which will modify our `SimpleSCM` in place:

```
julia> fit!(s)

Synthetic Control Model

Treatment panel:
Balanced Panel - single treated unit, continuous treatment
    Treated unit: United Kingdom
    Number of untreated units: 22
    First treatment period: 2016-07-01
    Number of pretreatment periods: 30
    Number of treatment periods: 9

        Model is fitted
        Impact estimates: [-0.54, -0.31, -0.206, -0.732, -1.241, -1.482, -1.818, -2.327, -1.994]
```

The reported impact estimates are the difference between observed outcome variable
and estimated outcome in the absence of treatment - a negative value therefore means
the treatment is expected to have reduced the outcome variable compared to the
counterfactual.

The package also defines a [plot recipe](https://github.com/JuliaPlots/RecipesBase.jl)
which allows to visualise the estimated impact:

```
julia> using Plots

julia> plot(s_model)
```
![Sample output](synthcontrol.png)

### Fitting a `SyntheticDiD` model

The package also implements a the synthetic differences-in-differences estimator of Arkhangelsky et
al. (2021) with the `SyntheticDiD` type. An example using data on German reunification:

```
julia> bp = load_germany_panel()
Balanced Panel - single treated unit, continuous treatment
    Treated unit: West Germany
    Number of untreated units: 16
    First treatment period: 1990
    Number of pretreatment periods: 30
    Number of treatment periods: 14
```

Here we are using the `load_*_panel()` family of functions rather than the `load_*()` family of
functions used above - when using the `panel` version of the `load` functions, a `BalancedPanel`
object is returned which obviates the need for creating this from the raw data. 

Fitting the model:

```
julia> sdid_model = SyntheticDiD(bp);

julia> fit!(sdid_model);
```

`show` and `plot` support is not yet implemented for `SyntheticDiD` models, but the point estimate
and standard error of the causal effect can be accessed directly from the model object:

```
julia> sdid_model.τ̂
1-element Vector{Float64}:
 -3121.95256130856

julia> sdid_model.se_τ̂
1-element Vector{Float64}:
 327.9157408571397
 ```
