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

The package also implements a the synthetic differences-in-differences estimator of [Arkhangelsky et
al. (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159) with the `SyntheticDiD` type.
An example using data on California's ban of tobacco advertising:

```
julia> sp = load_smoking_panel()
Balanced Panel - single treated unit, continuous treatment
    Treated unit: 3
    Number of untreated units: 38
    First treatment period: 1989
    Number of pretreatment periods: 19
    Number of treatment periods: 12
```

Here we are using the `load_*_panel()` family of functions rather than the `load_*()` family of
functions used above - when using the `panel` version of the `load` functions, a `BalancedPanel`
object is returned which obviates the need for creating this from the raw data. 

Fitting the model:

```
julia> sdid_model = SyntheticDiD(sp)
Synthetic Difference-in-Differences Model

Model is not fitted
```

As before, the `fit!` function is used to fit the model:

```
julia> fit!(sdid_model)
Synthetic Difference-in-Differences Model
        Model is fitted
        Impact estimate: -15.604
```

The model estimate can also be accessed as `sdid_model.τ̂`, and the standard error as
`sdid_model.se_τ̂`. The only algorithm for estimation of standard errors currently implemented is
the placebo algorithm, in which the estimator is sequentially applied to each control unit. By
default, standard errors are not estimated, the `se` keyword can be used to do so:

```
julia> fit!(sdid_model; se = :placebo)
Synthetic Difference-in-Differences Model
        Model is fitted
        Impact estimate: -15.604
                          (9.31)
```

### Fitting a `MC-NNM` model (experimental)

The package includes an experimental implementation of the Matrix Completion with Nuclear Norm
Minimization (MC-NNM) estimator ([Athey et al.,
2021](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1891924)). Due to its experimental
nature it is currently not exported and has to be accessed through the internal `fect_default`
function:

```
julia> SynthControl.fect_default(sp)
[ Info: Cross-validating...
        λ_norm = 1.0; MSPE = 184.11
        λ_norm = 0.4217; MSPE = 97.12
        λ_norm = 0.17783; MSPE = 60.18
        λ_norm = 0.07499; MSPE = 55.38
        λ_norm = 0.03162; MSPE = 83.68
        λ_norm = 0.01334; MSPE = 133.34
        λ_norm = 0.00562; MSPE = 176.03
        λ_norm = 0.00237; MSPE = 180.63
        λ_norm = 0.001; MSPE = 182.63
        λ_norm = 0.0; MSPE = 184.11

        λ_norm* = 0.07499
(τ̂ᵃᵗᵗ = -25.90741309106447, eff = [-3.777787420126515 -2.1951292344078013 … 2.8272388331494795 1.7451236037973956; -2.225234537444024 -0.6429018374729338 … -0.5174054676892013 -0.35710884584517544; … ; 2.3972577740355073 1.5448822989699096 … -1.6290392876871636 2.093390785866191; 1.1496735178409807 -1.3143906741400713 … -1.1176962789642886 -4.041039492499593])
```