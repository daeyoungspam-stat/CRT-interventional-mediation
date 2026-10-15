# crtmedsim

Reproducible simulations for causal mediation in cluster-randomized trials
with post-treatment confounding.

## Installation

From the repository root:

```sh
Rscript scripts/install_dependencies.R
R CMD INSTALL .
```

R 4.5 or later is required. The archived reference execution used R 4.5.1.
Install dependencies in a package library dedicated to the active R major/minor
version.
The package metadata uses a project-level placeholder contact. Release metadata
and public licensing terms can be finalized before publication.

## First run

```r
library(crtmedsim)

config <- simulation_config("smoke")
create_simulation("runs/smoke", config)
run_simulation("runs/smoke", cores = 1L)
bundle <- collect_results("runs/smoke")
results <- read_results(bundle)
head(results$summary)
```

The smoke preset runs one replication at K = 50. It retains the full fitting
settings, all scenarios, all methods, and both target scales. It checks execution;
Monte Carlo performance summaries require repeated replications.

## Full study

```r
create_simulation("runs/full", simulation_config("full"))
run_simulation("runs/full", cores = 10L)
bundle <- collect_results("runs/full")
```

The full preset has 2000 replications at each of K = 50, 75, 100, and 150.
Unix supports parallel replication workers. Windows runs with `cores = 1L`.
Repeating the same execution command resumes from completed matching files.
Workspace configuration is fixed at creation; a different configuration belongs
in a separate workspace.

For a selected execution range:

```r
run_simulation("runs/full", studies = "lowdim", K_values = 100L,
               rep_start = 1L, rep_end = 10L, cores = 10L)
```

Aggregation runs after all configured replication files from both studies are
available. Completion checks validate expected files, code and configuration
signatures, and duplicate result keys.

## Study contents

| Component | Setting |
|---|---|
| Targets | IDE, IIME, ISME |
| Scales | Cluster and individual |
| Cross-fitting | 5 treatment-stratified cluster folds |
| Canonical methods | Stabilized parametric EIF, stabilized ML EIF, parametric g-computation |
| Canonical scenarios | All 16 nuisance-status combinations |
| Low-dimensional references | Gaussian; residual bootstrap |
| Low-dimensional classifiers | glmnet; random forest; Super Learner |
| Low-dimensional scenarios | LD-AC, LD-T2 through LD-T6, LD-U2 through LD-U6 |
| Super Learner library | `c("SL.glm", "SL.ranger")` |
| Canonical ML integration | `max(80, K)` Gaussian draws per cluster |
| Reference draws | 10 per training cluster |
| G-computation intervals | 200 cluster bootstrap resamples |
| Classifier probability bounds | [0.025, 0.975] |
| Density and odds ratio bounds | [0.05, 20] |

The reference-classifier combinations give six low-dimensional methods.
The Gaussian reference generates fresh normal residuals. The residual-bootstrap
reference resamples centered, treatment-specific training residuals. The latter
is an empirical simulation implementation; its behavior is assessed through the
reported Monte Carlo results.

## Output

```text
runs/full/
  simulation_config.rds
  engine/
  runtime/
  logs/
  results/
    parts/canonical/K50/rep_0001.rds
    parts/lowdim/K50/rep_0001.rds
    annals_simulation_bundle.rds
    report/
      simulation_summary.csv
      canonical_scenarios.csv
      gcomp_scenarios.csv
      lowdim_scenarios.csv
      runtime_summary.csv
      truncation_summary.csv
      replication_warnings.csv
      comprehensive_simulation_tables.tex
      scenario_guarantees.tex
      figures/
```

The result bundle contains estimates, standard errors, intervals, target values,
scenario indicators, summaries, and diagnostics. Summary columns include bias,
empirical standard deviation, mean standard error, RMSE, coverage, and mean
confidence-interval width. Scenario expectation fields reproduce the archived
design's theoretical classifications and are distinct from empirical performance.

Existing compatible bundles can regenerate the archived report:

```r
render_report("results/annals_simulation_bundle.rds", "results/report")
```

## Reproduction checks

```r
validate_engine("runs/implementation-checks.log")
```

```sh
R CMD build .
R CMD check --no-manual crtmedsim_0.1.0.tar.gz
```

Set `CRTMEDSIM_INTEGRATION=true` when running `R CMD check` to also execute
the single-core, complete-method smoke study and check result aggregation,
report files, and resume behavior. This requires the fitting and report packages.

The engine files are matched to the final archived run by checksums. The package
tests check those hashes, configuration defaults, deterministic data generation,
target values, and workspace protection. `validate_engine()` executes additional
numerical implementation checks. Package and R versions affect fitted learners
and numerical results; each run records fitting-package versions. The recorded
reference bundle identifies R 4.5.1 and aggregation-time package versions.

See [design and implementation](docs/DESIGN.md),
[verification](docs/VERIFICATION.md), [server execution](docs/SLURM.md),
and [repository upload](docs/UPLOAD_KR.md).

