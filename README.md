# SynthControl.jl
Julia package for synthetic control methods, implementing the technique described in [Abadie et al. (2010)](https://economics.mit.edu/files/11859)

The package is an early draft and rather proof of concept.

## Installation

The package is currently unregistered, installation therefore works from directly
from the repo:

```
pkg> add "https://github.com/nilshg/SynthControl.jl"
```

## Usage

The package includes example data borrowed from the [CER Brexit study](
  https://www.cer.eu/insights/cost-brexit-june-2018):

```
julia> using SynthControl

julia> df = load_brexit()
897×4 DataFrames.DataFrame
│ Row │ country        │ qdate  │ realgdp │ dateid │
│     │ String         │ String │ Float64 │ Int64  │
├─────┼────────────────┼────────┼─────────┼────────┤
│ 1   │ Australia      │ 2009q1 │ 1.04    │ 57     │
│ 2   │ Austria        │ 2009q1 │ -1.53   │ 57     │
⋮
│ 895 │ Switzerland    │ 2018q3 │ 14.35   │ 95     │
│ 896 │ United Kingdom │ 2018q3 │ 15.72   │ 95     │
│ 897 │ United States  │ 2018q3 │ 19.32   │ 95     │
```

The package defines a `SynthControlModel` type, instances of which can be constructed
from a `DataFrame`, specifying the column holding the outcome variable of interest,
the time and group (ID) dimension, the treatment start period and the treated
group/observation. Currently, only one treatment unit can be specified.

The example data set includes quarterly GDP for a number of OECD countries, and
we are interested in estimating the impact of the Brexit vote in Q2 2016 on GDP
in the UK:

```
julia> s_model = SynthControlModel(df, :country, :dateid, :realgdp, 86, "United Kingdom")
Synthetic Control Model
Outcome variable: realgdp
Time dimension: dateid with 39 unique values
Treatment period: 86
ID variable: country with 23 unique values
Treatment ID: United Kingdom
Model is not fitted
```

The output indicates that the model is not fitted, that is we have at this stage
only defined the basic model structure. We can fit the model using the `fit!`
function, which will modify our `SynthControlModel` in place:

```
julia> fit!(s_model)
Synthetic Control Model

Outcome variable: realgdp
Time dimension: dateid with 39 unique values
Treatment period: 86
ID variable: country with 23 unique values
Treatment ID: United Kingdom
Model is fitted
Impact estimates: [-0.911, -1.053, -0.86, -0.672, -0.981, -1.509, -1.772, -2.165, -2.652, -2.25]
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
[Sample output](readme.png)

To do:
* Inference
* Matching on specific covariates
* Multiple treatment units
* Documentation
* Tests
* Consider dropping DataFrames dependency in favour of Tables
