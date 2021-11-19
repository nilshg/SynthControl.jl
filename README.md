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

The package includes example data borrowed from the [CER Brexit study](
  https://www.cer.eu/insights/cost-brexit-june-2018):

```
julia> using SynthControl

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

The package defines a `SynthControlModel` type, instances of which can be constructed
from a `TreatmentPanel` object from the package
[`TreatmentPanels`](https://github.com/nilshg/TreatmentPanels.jl). The `TreatmentPanel` is
constructed from a `DataFrame` and a specification of treatment assignment.

The example data set includes quarterly GDP for a number of OECD countries, and
we are interested in estimating the impact of the Brexit vote in Q2 2016 on GDP
in the UK:

```
julia> bp = BalancedPanel(df, "United Kingdom" => Date(2016, 7, 1); id_var = :country, t_var = :quarter, outcome_var = :realgdp)
Balanced Panel - single unit, single continuous treatment
    Treated unit: ["United Kingdom"]
    Number of untreated units: 22
    First treatment period: [Date("2016-07-01")]
    Number of pretreatment periods: [30]
    Number of treatment periods: [9]


julia> s = SynthControlModel(bp)

Synthetic Control Model

Treatment panel:
Balanced Panel - single unit, single continuous treatment
    Treated unit: ["United Kingdom"]
    Number of untreated units: 22
    First treatment period: [Date("2016-07-01")]
    Number of pretreatment periods: [30]
    Number of treatment periods: [9]


Model is not fitted
```

The output indicates that the model is not fitted, that is we have at this stage
only defined the basic model structure. We can fit the model using the `fit!`
function, which will modify our `SynthControlModel` in place:

```
julia> fit!(s)

Synthetic Control Model

Treatment panel:
Balanced Panel - single unit, single continuous treatment
    Treated unit: ["United Kingdom"]
    Number of untreated units: 22
    First treatment period: [Date("2016-07-01")]
    Number of pretreatment periods: [30]
    Number of treatment periods: [9]

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
