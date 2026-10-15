#!/usr/bin/env Rscript

source("annals_core.R")
source("annals_canonical.R")
source("annals_lowdim.R")

ap_assert(requireNamespace("data.table", quietly = TRUE),
          "data.table is required for post-processing.")

args <- ap_parse_args(commandArgs(trailingOnly = TRUE))
config <- ap_default_config()
config$reps <- as.integer(ap_arg(args, "reps", config$reps))
config$K_values <- as.integer(strsplit(ap_arg(
  args, "K-values", paste(config$K_values, collapse = ",")
), ",", fixed = TRUE)[[1L]])
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

output_dir <- ap_arg(args, "outdir", "final_results")
report_dir <- file.path(output_dir, "report")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
signature <- ap_config_signature(config)
code_signature <- ap_estimator_code_signature()

expected <- unlist(lapply(c("canonical", "lowdim"), function(study) {
  unlist(lapply(config$K_values, function(K) {
    vapply(seq_len(config$reps), function(rep_id) {
      ap_result_path(output_dir, study, K, rep_id)
    }, "")
  }), use.names = FALSE)
}), use.names = FALSE)
missing <- expected[!file.exists(expected)]
if (length(missing)) {
  writeLines(missing, file.path(report_dir, "missing_result_files.txt"))
  ap_stop(length(missing), " expected replication files are missing. See missing_result_files.txt.")
}

cores <- as.integer(ap_arg(args, "cores", 1L))
read_one <- function(path) {
  object <- readRDS(path)
  ap_assert(isTRUE(object$complete), paste("Incomplete result:", path))
  ap_assert(identical(object$config_signature, signature),
            paste("Config signature mismatch:", path))
  ap_assert(identical(object$code_signature, code_signature),
            paste("Estimator-code signature mismatch:", path))
  rows <- data.table::as.data.table(object$rows)
  rows[, elapsed_seconds := object$elapsed_seconds]
  diagnostics <- data.table::as.data.table(object$diagnostics)
  if (nrow(diagnostics)) {
    diagnostics <- diagnostics[, .(
      n = sum(n), truncated = sum(truncated),
      truncated_fraction = sum(truncated) / sum(n)
    ), by = .(study, K, rep, method, scenario, component, quantity, detail)]
  }
  warning_rows <- if (length(object$warnings)) {
    data.table::data.table(
      study = object$study, K = object$K, rep = object$rep,
      warning = object$warnings
    )
  } else data.table::data.table()
  list(rows = rows, diagnostics = diagnostics, warnings = warning_rows)
}
if (cores > 1L && .Platform$OS.type == "unix") {
  pieces <- parallel::mclapply(
    expected, read_one, mc.cores = cores, mc.preschedule = TRUE
  )
} else {
  pieces <- lapply(expected, read_one)
}
raw <- data.table::rbindlist(
  lapply(pieces, `[[`, "rows"), use.names = TRUE, fill = TRUE
)
diagnostics <- data.table::rbindlist(
  lapply(pieces, `[[`, "diagnostics"), use.names = TRUE, fill = TRUE
)
replication_warnings <- data.table::rbindlist(
  lapply(pieces, `[[`, "warnings"), use.names = TRUE, fill = TRUE
)
if (!ncol(replication_warnings)) {
  replication_warnings <- data.table::data.table(
    study = character(), K = integer(), rep = integer(), warning = character()
  )
}

key_columns <- c("study", "K", "rep", "scenario", "method", "estimand")
duplicate_keys <- raw[, .N, by = key_columns][N != 1L]
ap_assert(nrow(duplicate_keys) == 0L, "Duplicate result keys were detected.")
ap_assert(all(is.finite(raw$point)), "Non-finite point estimates were detected.")
ap_assert(all(is.finite(raw$true)), "Non-finite truth values were detected.")

summary <- raw[, .(
  n_rep = .N,
  truth = unique(true)[1L],
  mean_estimate = mean(point),
  bias = mean(error),
  bias_mcse = sd(error) / sqrt(.N),
  empirical_sd = sd(point),
  mean_se = mean(se),
  se_to_sd = mean(se) / sd(point),
  rmse = sqrt(mean(error^2)),
  coverage = mean(cover),
  coverage_mcse = sqrt(mean(cover) * (1 - mean(cover)) / .N),
  mean_ci_width = mean(ci_width),
  point_consistency_expected = unique(point_consistency_expected)[1L],
  asymptotic_linearity_expected = unique(asymptotic_linearity_expected)[1L],
  ci_validity_expected = unique(ci_validity_expected)[1L],
  ci_method = unique(ci_method)[1L],
  mean_elapsed_seconds = mean(elapsed_seconds)
), by = .(study, K, scenario, method, estimand, effect, scale)]

method_labels <- c(
  canonical_eif_par_stabilized = "Canonical EIF: parametric",
  canonical_eif_ml_stabilized = "Canonical EIF: machine learning",
  gcomp_parametric = "Parametric g-computation",
  ld_gaussian_glmnet_stabilized = "LD Gaussian + glmnet",
  ld_gaussian_ranger_stabilized = "LD Gaussian + random forest",
  ld_gaussian_sl_stabilized = "LD Gaussian + Super Learner",
  ld_resid_boot_glmnet_stabilized = "LD residual bootstrap + glmnet",
  ld_resid_boot_ranger_stabilized = "LD residual bootstrap + random forest",
  ld_resid_boot_sl_stabilized = "LD residual bootstrap + Super Learner"
)
summary[, method_label := unname(method_labels[method])]
raw[, method_label := unname(method_labels[method])]

canonical_theory <- data.table::rbindlist(lapply(
  ap_canonical_scenarios(), function(x) data.table::data.table(
    scenario = x$id, eta = x$eta, phi = x$phi,
    kappa = x$kappa, xi = x$xi,
    IDE_point = x$point_consistency_expected[["IDE"]],
    IIME_point = x$point_consistency_expected[["IIME"]],
    ISME_point = x$point_consistency_expected[["ISME"]],
    IDE_AL = x$asymptotic_linearity_expected[["IDE"]],
    IIME_AL = x$asymptotic_linearity_expected[["IIME"]],
    ISME_AL = x$asymptotic_linearity_expected[["ISME"]],
    IDE_CI = x$ci_validity_expected[["IDE"]],
    IIME_CI = x$ci_validity_expected[["IIME"]],
    ISME_CI = x$ci_validity_expected[["ISME"]]
  )
))
gcomp_theory <- data.table::rbindlist(lapply(
  ap_canonical_scenarios(), function(x) data.table::data.table(
    scenario = x$id, eta = x$eta, phi = x$phi,
    kappa = x$kappa, xi = x$xi,
    IDE_point = x$gcomp_point_expected[["IDE"]],
    IIME_point = x$gcomp_point_expected[["IIME"]],
    ISME_point = x$gcomp_point_expected[["ISME"]],
    IDE_AL = x$gcomp_point_expected[["IDE"]],
    IIME_AL = x$gcomp_point_expected[["IIME"]],
    ISME_AL = x$gcomp_point_expected[["ISME"]],
    IDE_CI = x$gcomp_point_expected[["IDE"]],
    IIME_CI = x$gcomp_point_expected[["IIME"]],
    ISME_CI = x$gcomp_point_expected[["ISME"]]
  )
))
lowdim_theory <- data.table::rbindlist(lapply(
  ap_lowdim_scenarios(), function(x) data.table::data.table(
    scenario = x$id,
    theta_eta = x$theta$eta, d = x$theta$d, s = x$theta$s,
    tau_eta = x$tau$eta, e = x$tau$e, s_minus = x$tau$s,
    IDE_point = x$point_consistency_expected[["IDE"]],
    IIME_point = x$point_consistency_expected[["IIME"]],
    ISME_point = x$point_consistency_expected[["ISME"]],
    IDE_AL = x$asymptotic_linearity_expected[["IDE"]],
    IIME_AL = x$asymptotic_linearity_expected[["IIME"]],
    ISME_AL = x$asymptotic_linearity_expected[["ISME"]],
    IDE_CI = x$ci_validity_expected[["IDE"]],
    IIME_CI = x$ci_validity_expected[["IIME"]],
    ISME_CI = x$ci_validity_expected[["ISME"]]
  )
))

runtime_by_rep <- unique(raw[, .(study, K, rep, elapsed_seconds)])
runtime <- runtime_by_rep[, .(
  replications = data.table::uniqueN(rep),
  mean_elapsed_seconds = mean(elapsed_seconds),
  median_elapsed_seconds = median(elapsed_seconds),
  total_worker_hours = sum(elapsed_seconds) / 3600
), by = .(study, K)]

truncation_summary <- diagnostics[, .(
  n = sum(n), truncated = sum(truncated),
  truncated_fraction = sum(truncated) / sum(n)
), by = .(study, K, method, scenario, component, quantity, detail)]

data.table::fwrite(summary, file.path(report_dir, "simulation_summary.csv"))
data.table::fwrite(canonical_theory, file.path(report_dir, "canonical_scenarios.csv"))
data.table::fwrite(gcomp_theory, file.path(report_dir, "gcomp_scenarios.csv"))
data.table::fwrite(lowdim_theory, file.path(report_dir, "lowdim_scenarios.csv"))
data.table::fwrite(runtime, file.path(report_dir, "runtime_summary.csv"))
data.table::fwrite(
  truncation_summary, file.path(report_dir, "truncation_summary.csv")
)
data.table::fwrite(
  replication_warnings, file.path(report_dir, "replication_warnings.csv")
)

code_files <- sort(list.files(
  pattern = "\\.(R|slurm|sh)$", full.names = TRUE
))
code_md5 <- tools::md5sum(code_files)
require_certification <- identical(
  tolower(ap_arg(args, "require-certification", "true")), "true"
)
certification_path <- ap_arg(
  args, "certification",
  if (require_certification) file.path("validation", "design_certification.rds")
  else ""
)
if (require_certification) {
  ap_assert(file.exists(certification_path),
            paste("Required design certification is missing:", certification_path))
}
certification <- if (nzchar(certification_path) && file.exists(certification_path)) {
  readRDS(certification_path)
} else NULL
if (!is.null(certification)) {
  ap_assert(identical(certification$config_signature, signature),
            "Design certification and simulation config signatures differ.")
  ap_assert(identical(certification$code_signature, code_signature),
            "Design certification and estimator-code signatures differ.")
}

bundle <- list(
  bundle_version = "annals-simulation-bundle-v1",
  complete = TRUE,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  config = config,
  config_signature = signature,
  code_signature = code_signature,
  truth = ap_truth(config),
  raw = raw,
  summary = summary,
  canonical_scenarios = canonical_theory,
  gcomp_scenarios = gcomp_theory,
  lowdim_scenarios = lowdim_theory,
  runtime = runtime,
  truncation_diagnostics_by_rep = diagnostics,
  truncation_summary = truncation_summary,
  replication_warnings = replication_warnings,
  design_certification = certification,
  code_md5 = code_md5,
  session_info = capture.output(sessionInfo())
)
bundle_path <- file.path(output_dir, "annals_simulation_bundle.rds")
temporary <- paste0(bundle_path, ".tmp-", Sys.getpid())
saveRDS(bundle, temporary, compress = "gzip")
ap_assert(file.rename(temporary, bundle_path), "Could not finalize result bundle.")

cat("Bundle written to ", bundle_path, "\n", sep = "")
cat("Rows: ", format(nrow(raw), big.mark = ","), "\n", sep = "")
