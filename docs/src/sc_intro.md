# Introduction to Synthetic Control Models

Synthetic Control Models are a class of causal inference models first introduced in [Abadie and
Gardeazabal (2003)](https://www.aeaweb.org/articles?id=10.1257/000282803321455188)[^1]. They infer
causal effects of treatments by creating a "synthetic" comparator unit as a weighted average of the
untreated observations to impute the missing potential outcomes for the treated unit in the
post-treatment period.

### Basic setup

The notation throughout the documentation (and, to the extent possible, the source code implementing
the various methods in this package) broadly follows [Arkhangelsky and Imbens
(2023)](https://arxiv.org/abs/2311.15458)[^2]. We observe $i = 1, 2, ..., N$ units for $t = 1, 2, ..., T$ periods. Outcomes are recorded in the $(N \times T)$ matrix $Y$:

$$Y = \begin{pmatrix}
    Y_{1,1} & Y_{1, 2} & \cdots & Y_{1, T} \\
    Y_{2,1} & \ddots   &        & \vdots   \\
    \vdots  &          &  \ddots& \vdots   \\
    Y_{N,1} & \cdots   & \cdots       & Y_{N, T}   
    \end{pmatrix}_{(N \times T)}$$

In every period, each unit can either be treated or untreated, with treatment status denoted by a
binary indicator $W_{it} \in \{0,1\}$. A number $N_{tr}$ units in the population are then exposed to
a treatment, leaving $N_{co}$ untreated units as poential comparator or control units. In the
simplest case, only one unit is treated at a point in time $T_0$, so that the observation window can
be divided into two parts: $T_{pre} \equiv 1, 2, ..., T_0$ and $T_{post} \equiv T_0 + 1, T_0 + 2,
..., T$. 

In potential outcome notation, each unit has two potential outcomes at each point in time,
$Y_{it}(0)$ if untreated and $Y_{it}(1)$ if treated. The causal effect of treatment is defined as
the difference between these two outcomes:

$$\tau_{it} \equiv Y_{it}(1) - Y_{it}(0)$$

As by construction only one of these outcomes is actually observable, the aim of causal inference is
to impute the missing outcome; in the case of treated units we are therefore looking to estimate the
outcome that would have obtained had the unit not been treated. 

### The synthetic control 

The synthetic control approach to solving the problem of imputing an unobserved potential outcome
was inspired by the comparative case studies often used in social sciences. In these comparative
case studies, the evolution of an outcome over time in a treated unit of interest (e.g. a country)
is compared to that of another treated unit that is somehow similar to the treated unit along some
relevant dimensions - essentially a form of matchin estimator. 

In practice however there often is not a single comparator unit which adequately approximates the
treated unit on the characteristica of interest, casting doubt on the validity of a simple
unit-to-unit comparison. Synthetic control methods address this issue by creating a synthetic
control unit which can be compared to the treated unit in such a way that its characteristics
closely resemble those of the treated unit.

Formally, let $X_{i, p}$ denote an $(N \times p)$ matrix of $p$ covariates (which can contain
pre-treatment outcomes) measured for each of the $N$ observed units, and let $\omega$ denote an
$(N_{co} \times 1)$ vector of weights for each untreated unit. Synthetic control methods find the
weights $\omega$ by solving an optimization problem of the form 

$$\omega^* = \arg \min_{\omega}|| X_{tr, \cdot} - X_{co, \cdot}'\omega||_2$$

Generally, this optimization problem is constrained in some way, e.g. by requiring the weights
$\omega$ to lie between 0 and 1 for each unit and to sum to 1 across all control units, with
different synthetic control estimators differing in their exact specification of this optimization
problem and its constraints. 

With the optimal weights $\omega^*$ in hand, imputation of the missing potential outcomes for the
treated unit is done by simply calculating the average outcome across control units, weighted by
$\omega^*$:

$$\hat{\tau}_{sc} = Y_{tr, T_{post}} - Y_{co, T_{post}}'\omega^*$$

[^1]: Abadie, A., and Gardeazabal, J. (2003): *The Economic Cost of Conflict: A Case Study of the
    Basque Country*, American Economic Review, Vol. 93(1), Pp. 113-132
[^2]: Arkhangelsky, D., and Imbens, G. (2023): *Causal Models for Longitudinal and Panel Data: A
    Survey*, arXiv:2311.15458