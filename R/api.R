simulation_config <- function(preset = c("full", "smoke")) {
  preset <- match.arg(preset)
  config <- .core()$ap_default_config()
  if (preset == "smoke") {
    config$reps <- 1L
    config$K_values <- 50L
  }
  config
}

simulate_trial <- function(K = 50L, rep_id = 1L,
                           config = simulation_config()) {
  config <- .validate_config(config)
  K <- .integer_vector(K, "K", 20L, TRUE)
  rep_id <- .integer_vector(rep_id, "rep_id", scalar = TRUE)
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  .core()$ap_generate_data(K, config, rep_id, include_y = TRUE)
}

simulation_truth <- function(config = simulation_config()) {
  .core()$ap_truth(.validate_config(config))
}

create_simulation <- function(path, config = simulation_config("smoke")) {
  config <- .validate_config(config)
  if (dir.exists(path) && length(list.files(path, all.files = TRUE,
                                           no.. = TRUE)) &&
      !file.exists(file.path(path, "simulation_config.rds"))) {
    stop("Choose an empty directory or an existing simulation workspace.", call. = FALSE)
  }
  path <- .absolute_dir(path, create = TRUE)
  config_path <- file.path(path, "simulation_config.rds")
  if (file.exists(config_path) && !identical(readRDS(config_path), config)) {
    stop("Workspace configuration differs; choose a separate directory.", call. = FALSE)
  }
  .copy_engine(file.path(path, "engine"))
  if (!file.exists(config_path)) saveRDS(config, config_path)
  for (name in c("results", "logs", "runtime")) {
    dir.create(file.path(path, name), showWarnings = FALSE)
  }
  invisible(path)
}

run_simulation <- function(path, studies = c("canonical", "lowdim"),
                           K_values = NULL, rep_start = 1L, rep_end = NULL,
                           cores = 1L, overwrite = FALSE) {
  workspace <- .load_workspace(path)
  path <- workspace$path
  config <- workspace$config
  if (!is.character(studies) || !length(studies) || anyNA(studies) ||
      anyDuplicated(studies) || !all(studies %in% c("canonical", "lowdim"))) {
    stop("studies must select canonical and/or lowdim.", call. = FALSE)
  }
  if (is.null(K_values)) K_values <- config$K_values
  K_values <- .integer_vector(K_values, "K_values", 20L)
  if (!all(K_values %in% config$K_values)) stop("K_values exceeds the workspace design.")
  if (is.null(rep_end)) rep_end <- config$reps
  rep_start <- .integer_vector(rep_start, "rep_start", scalar = TRUE)
  rep_end <- .integer_vector(rep_end, "rep_end", scalar = TRUE)
  cores <- .integer_vector(cores, "cores", scalar = TRUE)
  if (rep_start > rep_end || rep_end > config$reps) stop("Invalid replication range.")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  if (.Platform$OS.type == "windows" && cores > 1L) {
    stop("Use cores = 1 on Windows; multicore execution requires Unix.")
  }
  packages <- c("SuperLearner", "glmnet", "ranger")
  .dependencies(packages)
  versions <- data.frame(
    package = c("R", "crtmedsim", packages),
    version = c(as.character(getRversion()), as.character(utils::packageVersion("crtmedsim")),
                vapply(packages, function(x) as.character(utils::packageVersion(x)), ""))
  )
  for (study in studies) for (K in K_values) {
    tag <- sprintf("%s_K%d_reps%d-%d", study, K, rep_start, rep_end)
    log <- file.path(path, "logs", paste0(tag, ".log"))
    utils::write.csv(versions, file.path(path, "runtime", paste0(tag, "_versions.csv")),
                     row.names = FALSE)
    message(sprintf("%s: K=%d, replications %d-%d, cores=%d", study, K,
                    rep_start, rep_end, cores))
    .execute(file.path(path, "engine"), "annals_run.R", c(
      .config_args(config), "--study", study, "--K", K,
      "--rep-start", rep_start, "--rep-end", rep_end, "--cores", cores,
      "--outdir", file.path(path, "results"),
      "--overwrite", tolower(as.character(overwrite))
    ), log)
  }
  invisible(file.path(path, "results"))
}

collect_results <- function(path, figures = TRUE) {
  workspace <- .load_workspace(path)
  path <- workspace$path
  .dependencies("data.table")
  .execute(file.path(path, "engine"), "annals_postprocess.R", c(
    .config_args(workspace$config), "--outdir", file.path(path, "results"),
    "--cores", "1", "--require-certification", "false"
  ), file.path(path, "logs", "collect.log"))
  bundle <- file.path(path, "results", "annals_simulation_bundle.rds")
  if (isTRUE(figures)) render_report(bundle, file.path(path, "results", "report"))
  invisible(bundle)
}

render_report <- function(bundle, output_dir) {
  .dependencies(c("data.table", "ggplot2"))
  bundle <- normalizePath(bundle, winslash = "/", mustWork = TRUE)
  output_dir <- .absolute_dir(output_dir, create = TRUE)
  engine <- tempfile("crtmedsim-engine-")
  on.exit(unlink(engine, recursive = TRUE), add = TRUE)
  .copy_engine(engine)
  .execute(engine, "annals_make_report.R", c(
    "--bundle", bundle, "--output-dir", output_dir,
    "--require-certification", "false"
  ), file.path(output_dir, "render.log"))
  invisible(output_dir)
}

read_results <- function(bundle) {
  object <- readRDS(bundle)
  expected <- .core()$ap_estimator_code_signature(.engine_dir())
  if (!isTRUE(object$complete) || !identical(object$code_signature, expected)) {
    stop("A complete bundle with the matching engine signature is required.", call. = FALSE)
  }
  object
}

validate_engine <- function(log = tempfile("crtmedsim-validation-", fileext = ".log")) {
  .dependencies(c("SuperLearner", "glmnet", "ranger"))
  log <- file.path(.absolute_dir(dirname(log), create = TRUE), basename(log))
  engine <- tempfile("crtmedsim-engine-")
  on.exit(unlink(engine, recursive = TRUE), add = TRUE)
  .copy_engine(engine)
  .execute(engine, "validate_implementation.R", log = log)
  invisible(log)
}

