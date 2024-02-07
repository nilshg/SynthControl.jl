## Simple Synthetic Control Model

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