library(crtmedsim)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript scripts/reproduce.R smoke|full WORKSPACE [CORES]")
}
preset <- match.arg(args[[1L]], c("smoke", "full"))
cores <- if (length(args) == 3L) as.integer(args[[3L]]) else 1L
create_simulation(args[[2L]], simulation_config(preset))
run_simulation(args[[2L]], cores = cores)
collect_results(args[[2L]])

