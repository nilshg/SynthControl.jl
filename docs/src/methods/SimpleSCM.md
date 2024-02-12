# Simple Synthetic Control Model

This method captures the simple-most version of a Synthetic Control Model. It only relies on
outcome data ``Y_{it}`` and treatment indicators ``W_{it}`` to find unit-weights ``\omega_i`` by solving

$$\omega^* = \arg \min_{\omega} || Y_{tr, 1:T_0} - \omega Y_{co, 1:T_0} ||_2 \\
s.t. \sum_{i = 1}^{N_{co}} \omega_i = 1 \\
0 \leq \omega_i \leq 1 \ \forall \ i \in 1, ..., N_{co}$$

where $Y_{tr, \cdot}$ is a length $T_0$ vector of pre-treatment outcomes for the treated unit,
$Y_{co, 1:T_0}$ is a $(N_{co} \times T_0)$ matrix of pre-treatment outcomes for the $N_{co}$ control
units and $\omega$ is a length $N_{co}$ vector of weights for each control unit. Weights are
restricted to lie between zero and one for each control unit and to sum up to one, as in Abadie and
Gardeazabal (2003). 

The `SimpleSCM` method can therefore be seen as a simplified version of the Abadie and Gardeazabal
(2003) or Abadie, Diamond and Hainmüller (2011) Synthetic Control Model, with the following
differences:

* Only past outcomes are being used to find weights, no other covariate information is taken into
  account;
* There is no weighting of covariates (i.e., past outcomes) - in the original Abadie/Gardeazabal
  notation, the diagonal matrix $V$ which weigts covariates is the identity matrix. 

Despite these limitations, the results from this simple procedure can often be close to those of the
more fully featured `ADH2010` model, as shown in the example below. 

## Implementation in `SynthControl`

The simple SCM model can be estimated by passing a `TreatmentPanel` to the `SimpleSCM` constructor.
The following example estimates the effect of California's proposition 99 on cigarette sales, an
application taken from [Abadie, Diamond and Hainmüller (2010,
JASSA)](https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.ap08746)[^1].

First, we load the relevant data:

```julia
julia> using SynthControl

julia> smoking_panel = load_smoking_panel()
Balanced Panel - single treated unit, continuous treatment
    Treated unit: 3
    Number of untreated units: 38
    First treatment period: 1989
    Number of pretreatment periods: 19
    Number of treatment periods: 12
```

Then we create the `SimpleSCM` model and call `fit!` on the model instance:

```julia
julia> simple_scm = SimpleSCM(smoking_panel);

julia> fit!(simple_scm)
Synthetic Control Model

Treatment panel:
Balanced Panel - single treated unit, continuous treatment
    Treated unit: 3
    Number of untreated units: 38
    First treatment period: 1989
    Number of pretreatment periods: 19
    Number of treatment periods: 12

        Model is fitted
        Impact estimates: [-8.44, -9.207, -12.634, -13.729, -17.534, -22.049, -22.858, -23.997, -26.261, -23.338, -27.52, -26.597]
        ATT: -19.514
```

The average treatment effect on the treated (ATT) is simply the average of the imputed synthetic
control outcomes for the post-treatment periods:

```
julia> using Statistics

julia> mean(simple_scm.τ̂)
-19.51362976399461
```

This 

[^1]: Abadie, A., Diamond, A., and Hainmüller, J. (2010): *Synthetic Control Methods for Comparative
    Case Studies: Estimating the Effect of California’s Tobacco Control Program*, American Journal
    of Political Science, Vol. 59(2), Pp. 495-510