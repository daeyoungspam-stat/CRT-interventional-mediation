library(crtmedsim)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: array.R WORKSPACE STUDY REPLICATIONS_PER_TASK")
workspace <- args[[1L]]
study <- match.arg(args[[2L]], c("canonical", "lowdim"))
chunk <- as.integer(args[[3L]])
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
stopifnot(length(chunk) == 1L, is.finite(chunk), chunk >= 1L,
          is.finite(task), task >= 1L)
config <- readRDS(file.path(workspace, "simulation_config.rds"))
per_K <- ceiling(config$reps / chunk)
index <- (task - 1L) %/% per_K + 1L
if (index > length(config$K_values)) stop("Task ID exceeds the design.")
start <- ((task - 1L) %% per_K) * chunk + 1L
end <- min(config$reps, start + chunk - 1L)
run_simulation(workspace, studies = study, K_values = config$K_values[[index]],
               rep_start = start, rep_end = end, cores = cores)

