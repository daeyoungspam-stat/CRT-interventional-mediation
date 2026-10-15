\name{crtmedsim}
\alias{crtmedsim}
\alias{simulation_config}
\alias{simulate_trial}
\alias{simulation_truth}
\alias{create_simulation}
\alias{run_simulation}
\alias{collect_results}
\alias{render_report}
\alias{read_results}
\alias{validate_engine}
\title{Reproduce the Fixed Cluster-Randomized Mediation Simulation}
\description{
Configure and execute the archived simulation engine, assemble replication
results, and produce numerical tables and figures.
}
\usage{
simulation_config(preset = c("full", "smoke"))
simulate_trial(K = 50L, rep_id = 1L, config = simulation_config())
simulation_truth(config = simulation_config())
create_simulation(path, config = simulation_config("smoke"))
run_simulation(path, studies = c("canonical", "lowdim"),
  K_values = NULL, rep_start = 1L, rep_end = NULL,
  cores = 1L, overwrite = FALSE)
collect_results(path, figures = TRUE)
render_report(bundle, output_dir)
read_results(bundle)
validate_engine(log = tempfile("crtmedsim-validation-", fileext = ".log"))
}
\arguments{
\item{preset}{\code{"full"} selects 2000 replications at 50, 75, 100, and
150 clusters. \code{"smoke"} selects one replication at 50 clusters with the
same fitting and integration settings.}
\item{K}{Number of clusters in a simulated dataset.}
\item{rep_id}{Positive integer replication identifier controlling the data seed.}
\item{config}{Configuration list returned by \code{simulation_config()}.}
\item{path}{Simulation workspace directory. Create it with \code{create_simulation()}.}
\item{studies}{Selection from \code{"canonical"} and \code{"lowdim"}.}
\item{K_values}{Cluster counts to execute from the workspace design; \code{NULL}
selects all configured counts.}
\item{rep_start, rep_end}{Inclusive replication range; a \code{NULL} end selects
the configured replication count.}
\item{cores}{Replication workers. The default is one. Multicore execution uses
forked workers on Unix; Windows supports one worker.}
\item{overwrite}{Recompute completed replications when \code{TRUE}.}
\item{figures}{Generate figures and LaTeX tables after aggregation.}
\item{bundle}{Path to an engine-compatible result bundle in RDS format.}
\item{output_dir}{Directory for regenerated reports.}
\item{log}{Path for the implementation-check log.}
}
\details{
The engine includes stabilized canonical parametric and machine-learning EIF
estimators, parametric g-computation, and six low-dimensional reference-classifier
combinations. Targets comprise IDE, IIME, and ISME on cluster and individual scales.
The residual-bootstrap reference resamples centered training residuals and is
an empirical implementation evaluated by the simulation.

Runtime dependencies for fitting are \pkg{SuperLearner}, \pkg{glmnet}, and
\pkg{ranger}. Aggregation uses \pkg{data.table}; reports additionally use
\pkg{ggplot2}. Subprocesses inherit the active R library paths and limit numerical
library threads to one per worker. Replication seeds depend on identifiers.

Aggregation requires all configured replications of both studies. Saved
replications include configuration and estimator-code signatures. An existing
matching complete replication is reused when \code{overwrite = FALSE}.
Each workspace stores the engine, configuration, runtime package versions,
replication files, and logs.
}
\value{
\code{simulation_config()} returns a list. \code{simulate_trial()} returns a
long-format data frame containing generated observations and simulation-specific
conditional-mean columns. \code{simulation_truth()} returns a named numeric vector.
\code{read_results()} returns the saved bundle list. Other functions invisibly
return a directory, bundle path, or log path.
}
\examples{
config <- simulation_config("smoke")
trial <- simulate_trial(K = 50L, rep_id = 1L, config = config)
head(trial)
simulation_truth(config)
\dontrun{
create_simulation("runs/smoke", config)
run_simulation("runs/smoke", cores = 1L)
bundle_path <- collect_results("runs/smoke")
results <- read_results(bundle_path)
head(results$summary)
}
}

