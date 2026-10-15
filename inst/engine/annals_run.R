#!/usr/bin/env Rscript

source("annals_core.R")
source("annals_canonical.R")
source("annals_lowdim.R")

ap_required_packages <- function(study) {
  if (study == "canonical") {
    return(c("SuperLearner", "glmnet", "ranger"))
  }
  c("SuperLearner", "glmnet", "ranger")
}

ap_assert_runtime_packages <- function(study) {
  packages <- ap_required_packages(study)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    ap_stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

ap_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    ap_stop("Atomic rename failed for ", path)
  }
  invisible(path)
}

ap_run_one_rep <- function(study, K, rep_id, config, output_dir,
                           overwrite = FALSE) {
  path <- ap_result_path(output_dir, study, K, rep_id)
  signature <- ap_config_signature(config)
  code_signature <- ap_estimator_code_signature()
  if (file.exists(path) && !overwrite) {
    existing <- try(readRDS(path), silent = TRUE)
    if (!inherits(existing, "try-error") &&
        identical(existing$config_signature, signature) &&
        identical(existing$code_signature, code_signature) &&
        isTRUE(existing$complete)) {
      return(list(path = path, skipped = TRUE, elapsed = 0))
    }
  }
  started <- proc.time()[["elapsed"]]
  warning_messages <- character(0)
  output <- withCallingHandlers({
    data <- ap_generate_data(K, config, rep_id, include_y = TRUE)
    rows <- diagnostics <- list()
    if (study == "canonical") {
      if ("par" %in% config$canonical_estimators) {
        fitted <- ap_run_canonical_eif(
          data, K, rep_id, "par", config, stabilize = TRUE
        )
        rows$par <- fitted$rows
        diagnostics$par <- fitted$diagnostics
      }
      if ("ml" %in% config$canonical_estimators) {
        fitted <- ap_run_canonical_eif(
          data, K, rep_id, "ml", config, stabilize = TRUE
        )
        rows$ml <- fitted$rows
        diagnostics$ml <- fitted$diagnostics
      }
      if ("gcomp" %in% config$canonical_estimators) {
        rows$gcomp <- ap_run_gcomp(data, K, rep_id, config)
      }
    } else if (study == "lowdim") {
      fitted <- ap_run_lowdim_eif(
        data, K, rep_id, config, stabilize = TRUE
      )
      rows$lowdim <- fitted$rows
      diagnostics$lowdim <- fitted$diagnostics
    } else {
      ap_stop("study must be canonical or lowdim.")
    }
    list(
      rows = do.call(rbind, rows),
      diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else
        data.frame()
    )
  }, warning = function(warning) {
    warning_messages <<- c(warning_messages, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  elapsed <- proc.time()[["elapsed"]] - started
  expected_rows <- ap_expected_rows_per_rep(study, config)
  ap_assert(nrow(output$rows) == expected_rows, sprintf(
    "%s replication %d at K=%d produced %d rows; expected %d.",
    study, rep_id, K, nrow(output$rows), expected_rows
  ))
  result_key <- paste(
    output$rows$scenario, output$rows$method, output$rows$estimand,
    sep = "|"
  )
  ap_assert(!anyDuplicated(result_key),
            "Duplicate scenario-method-estimand rows were produced.")
  numeric_output <- c("point", "se", "ci_low", "ci_up", "true")
  ap_assert(all(vapply(output$rows[numeric_output], function(x) {
    all(is.finite(x))
  }, logical(1))), "A replication produced non-finite inferential output.")
  result <- list(
    complete = TRUE,
    study = study,
    K = K,
    rep = rep_id,
    rows = output$rows,
    diagnostics = output$diagnostics,
    elapsed_seconds = elapsed,
    config_signature = signature,
    code_signature = code_signature,
    seed = ap_hash_seed(config$base_seed, "dgp", K, rep_id),
    truth = ap_truth(config),
    warnings = unique(warning_messages),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  ap_atomic_save_rds(result, path)
  list(path = path, skipped = FALSE, elapsed = elapsed,
       warnings = unique(warning_messages))
}

args <- ap_parse_args(commandArgs(trailingOnly = TRUE))
study <- ap_arg(args, "study", Sys.getenv("STUDY", ""))
ap_assert(study %in% c("canonical", "lowdim"),
          "Use --study canonical or --study lowdim.")

config <- ap_default_config()
config$reps <- as.integer(ap_arg(args, "reps", config$reps))
config$K_values <- as.integer(strsplit(
  ap_arg(args, "K-values", paste(config$K_values, collapse = ",")), ",",
  fixed = TRUE
)[[1L]])
config$folds <- as.integer(ap_arg(args, "folds", config$folds))
config$mc_draws_ml <- as.integer(ap_arg(args, "mc-draws-ml", config$mc_draws_ml))
config$gcomp_boot <- as.integer(ap_arg(args, "gcomp-boot", config$gcomp_boot))
config$reference_draws <- as.integer(ap_arg(
  args, "reference-draws", config$reference_draws
))
config$base_seed <- as.integer(ap_arg(args, "base-seed", config$base_seed))
config$canonical_estimators <- strsplit(ap_arg(
  args, "canonical-estimators",
  paste(config$canonical_estimators, collapse = ",")
), ",", fixed = TRUE)[[1L]]
config$lambda_laws <- strsplit(ap_arg(
  args, "lambda-laws", paste(config$lambda_laws, collapse = ",")
), ",", fixed = TRUE)[[1L]]
config$classifier_names <- strsplit(ap_arg(
  args, "classifier-names", paste(config$classifier_names, collapse = ",")
), ",", fixed = TRUE)[[1L]]
ap_validate_config(config)
ap_assert_runtime_packages(study)

output_dir <- ap_arg(args, "outdir", "final_results")
cores <- as.integer(ap_arg(args, "cores", 1L))
task_reps <- as.integer(ap_arg(args, "task-reps", 10L))
overwrite <- identical(tolower(ap_arg(args, "overwrite", "false")), "true")

explicit_K <- ap_arg(args, "K", NULL)
explicit_start <- ap_arg(args, "rep-start", NULL)
explicit_end <- ap_arg(args, "rep-end", NULL)
if (!is.null(explicit_K)) {
  K <- as.integer(explicit_K)
  rep_start <- as.integer(if (is.null(explicit_start)) 1L else explicit_start)
  rep_end <- as.integer(if (is.null(explicit_end)) rep_start else explicit_end)
} else {
  task_id <- as.integer(ap_arg(
    args, "task-id", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")
  ))
  chunks_per_K <- ceiling(config$reps / task_reps)
  K_index <- (task_id - 1L) %/% chunks_per_K + 1L
  ap_assert(K_index <= length(config$K_values), "Array task id exceeds design.")
  chunk_index <- (task_id - 1L) %% chunks_per_K + 1L
  K <- config$K_values[[K_index]]
  rep_start <- (chunk_index - 1L) * task_reps + 1L
  rep_end <- min(config$reps, rep_start + task_reps - 1L)
}
ap_assert(K %in% config$K_values, "Requested K is not in config$K_values.")
ap_assert(rep_start >= 1L && rep_end <= config$reps && rep_start <= rep_end,
          "Invalid replication range.")

rep_ids <- seq.int(rep_start, rep_end)
cat(sprintf(
  "Study=%s K=%d reps=%d-%d cores=%d signature=%s\n",
  study, K, rep_start, rep_end, cores, ap_config_signature(config)
))

worker <- function(rep_id) {
  tryCatch(
    ap_run_one_rep(study, K, rep_id, config, output_dir, overwrite),
    error = function(error) structure(list(
      rep = rep_id, message = conditionMessage(error)
    ), class = "ap_rep_error")
  )
}

if (cores > 1L && .Platform$OS.type == "unix") {
  results <- parallel::mclapply(
    rep_ids, worker, mc.cores = min(cores, length(rep_ids)),
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  results <- lapply(rep_ids, worker)
}

failed <- vapply(results, inherits, logical(1), "ap_rep_error")
if (any(failed)) {
  for (failure in results[failed]) {
    message("Replication ", failure$rep, " failed: ", failure$message)
  }
  ap_stop(sum(failed), " replication(s) failed.")
}
elapsed <- sum(vapply(results, function(x) x$elapsed, 0.0))
skipped <- sum(vapply(results, function(x) x$skipped, FALSE))
warnings <- unique(unlist(lapply(results, function(x) x$warnings), use.names = FALSE))
if (length(warnings)) {
  cat("Unique warnings captured by this task:\n")
  for (warning in warnings) cat("- ", warning, "\n", sep = "")
}
cat(sprintf(
  "Task complete: %d replications, %d skipped, %.1f worker-seconds.\n",
  length(results), skipped, elapsed
))
