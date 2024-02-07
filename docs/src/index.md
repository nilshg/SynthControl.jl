# SynthControl.jl 

_A Julia package for estimating causal effects via synthetic control estimators_

## Package Overview

`SynthControl` aims to provide a wide coverage of synthetic control estimators, written in pure
Julia. The package is under development, with new estimators and options for statistical inference
added. Where possible, estimators are implemented following reference implementations by the
original authors.

The following table provides an overview of the planned scope of the package as well as current
implementation status 

| Estimator |  Point estimate  | Covariates |  Standard Errors  |  Reference implementation |
|--------------|:----:|:----:|:-----:|------|
| Simple SCM   |  🟢  |  🟥  |  🟥  | None |
| ADH2010      |  🟡  |  🟥  |  🟥  | None |
| SyntheticDiD |  🟢  |  🟥  |  🟡  | [synthdid (R)](https://github.com/synth-inference/synthdid) |
| PenalizedSCM |  🟥  |  🟥  |  🟥  | [pensynth (R)](https://github.com/jeremylhour/pensynth) |
| AugmentedSCM |  🟥  |  🟥  |  🟥  | [augsynth (R)](https://github.com/ebenmichael/augsynth) |  
| MC-NNM       |  🟡  |  🟥  |  🟥  | [fect (R)](https://github.com/xuyiqing/fect/)  |

Contributions to add methods to the package scope - or even better, full implementations! - are very
much welcome.

## Documentation outline

* Introduction to Synthetic Control estimation
* Package design
* Available estimators
    * SimpleSCM
    * ADH2010
    * SyntheticDiD
    * PenalizedSCM
    * AugmentedSCM
    * MC-NNM
* Examples

## API 

```@docs
SimpleSCM
SyntheticDiD
fit!
isfitted
```