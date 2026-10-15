# Design and implementation

## Archived engine

The engine version is `annals-v8-sl-glm-ranger-dgp9`, with estimator signature
`estimator-code-0717194216`. `inst/engine/checksums.csv` records the seven source
file hashes. Each installed-package API call verifies these hashes. The interface
copies the engine into a workspace and invokes the existing execution scripts.

## Data generation

Clusters are independent. Cluster size is discrete uniform on 10 through 20.
The generator draws cluster and individual baseline covariates, assigns treatment
with probability 0.5, and generates exchangeable Gaussian processes for the
post-treatment confounder and mediator. The mediator conditional mean depends on
the confounder through coefficient 0.10. The process correlation is 0.06.
The outcome conditional mean includes individual and peer mediator/confounder
terms and interactions. The configuration records every numerical coefficient;
the generator and feature definitions are in `annals_core.R`.

```r
library(crtmedsim)
simulation_config()$dgp
simulation_truth()
```

The target values are:

| Target | Cluster scale | Individual scale |
|---|---:|---:|
| IDE | 0.468017950368579 | 0.489481156033335 |
| IIME | 0.330550604999998 | 0.342964557000000 |
| ISME | 0.192445875000001 | 0.188840114999999 |

## Canonical estimation

`annals_canonical.R` fits outcome and process models on training clusters and
evaluates the estimators on held-out clusters. Parametric fits use the archived
regression features. ML conditional means use the two-component Super Learner
library. Gaussian residual laws and estimated exchangeable correlations supply
the process densities. ML partial means use common Gaussian integration draws.

The correct parametric marginal mediator law is derived from the full confounder
and conditional mediator fits. These full fits are retained when a scenario
selects a reduced confounder or conditional mediator fit elsewhere. The correct
ML marginal mediator mean is fitted directly. Reduced process models retain
treatment and cluster size; the reduced conditional mediator model also retains
the confounder. Reduced outcome regressions retain treatment-specific process
terms while dropping baseline terms and their interactions. The source feature
functions specify each term exactly.

G-computation integrates the parametric outcome regression under the corresponding
fitted Gaussian laws and forms cluster-percentile bootstrap intervals.

## Low-dimensional estimation

`annals_lowdim.R` constructs full and mediator-only augmented classifiers and
evaluates their fitted odds ratio. The reference mean regresses the mediator on
cluster size, the cluster covariate, and the individual-covariate cluster mean,
separately within each treatment arm. Gaussian references add fresh normal draws;
residual-bootstrap references add resampled centered training residuals.

The classifier choices are cross-validated glmnet, random forest, and Super
Learner. Auxiliary process fits and regressions, scenario-specific feature
restrictions, marginal mediator ratios, and the treatment classifiers are
implemented in the same source file. The Super Learner library is
`c("SL.glm", "SL.ranger")`. Its logistic GLM component starts with unpenalized
fitting and falls back to ridge under the archived separation or convergence
criteria; the penalty is the inverse training-cluster count.

The empirical residual-bootstrap reference is retained as the archived
simulation procedure. Its interpretation is empirical, alongside the Gaussian
reference results.

## Stabilization and intervals

The EIF implementations update the outcome regression using a fitted residual
correction for each component, then recompute the component contributions.
Classifier probabilities and ratios receive the fixed bounds in the configuration.
The contrast contributions determine cluster-level standard errors. For the
individual scale, the calculation centers by the estimated target times cluster
size and divides by the mean cluster size. The interval multiplier is
`sqrt(K / (K - 3)) * qt(0.975, K - 3)` at the default confidence level.

## Randomness and result integrity

Seeds are derived from the base seed, replication identifier, cluster count, and
component identifiers. Both studies share the replication-specific generated
dataset. Workers execute complete replications. Atomic replication-file writes,
configuration signatures, and source signatures support interrupted runs.

Runtime fitting-package versions are written per execution range. Exact numerical
matching additionally depends on R, package versions, and numerical libraries.

