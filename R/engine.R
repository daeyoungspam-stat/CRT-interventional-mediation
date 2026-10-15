.engine_dir <- function() {
  path <- system.file("engine", package = "crtmedsim", mustWork = TRUE)
  manifest <- utils::read.csv(file.path(path, "checksums.csv"),
                              stringsAsFactors = FALSE)
  actual <- unname(tools::md5sum(file.path(path, manifest$file)))
  if (!identical(actual, manifest$md5)) {
    stop("Installed engine checksum mismatch.", call. = FALSE)
  }
  path
}

.core <- function() {
  env <- new.env(parent = asNamespace("crtmedsim"))
  old <- options("stringsAsFactors")
  on.exit(options(old), add = TRUE)
  sys.source(file.path(.engine_dir(), "annals_core.R"), envir = env)
  env
}

.integer_vector <- function(x, name, lower = 1L, scalar = FALSE) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x < lower | x > .Machine$integer.max | x != floor(x)) ||
      (scalar && length(x) != 1L) || anyDuplicated(x)) {
    stop(name, " must contain distinct integers >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

.validate_config <- function(config) {
  defaults <- .core()$ap_default_config()
  if (!is.list(config) || !identical(names(config), names(defaults))) {
    stop("Start with simulation_config() and retain its named fields.", call. = FALSE)
  }
  editable <- c("reps", "K_values", "folds", "base_seed", "mc_draws_ml",
                "gcomp_boot", "reference_draws", "canonical_estimators",
                "lambda_laws", "classifier_names")
  for (name in setdiff(names(defaults), editable)) {
    if (!identical(config[[name]], defaults[[name]])) {
      stop("The reproduction engine fixes configuration field: ", name,
           call. = FALSE)
    }
  }
  for (name in c("reps", "folds", "base_seed", "mc_draws_ml", "gcomp_boot",
                 "reference_draws")) {
    lower <- if (name %in% c("folds", "gcomp_boot")) 2L else 1L
    config[[name]] <- .integer_vector(config[[name]], name, lower, TRUE)
  }
  config$K_values <- .integer_vector(config$K_values, "K_values", 20L)
  if (config$folds >= min(config$K_values)) {
    stop("folds must be smaller than every cluster count.", call. = FALSE)
  }
  for (name in c("canonical_estimators", "lambda_laws", "classifier_names")) {
    value <- config[[name]]
    if (!is.character(value) || !length(value) || anyNA(value) ||
        anyDuplicated(value) || !all(value %in% defaults[[name]])) {
      stop("Invalid method selection in ", name, ".", call. = FALSE)
    }
  }
  .core()$ap_validate_config(config)
  config
}

.absolute_dir <- function(path, create = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Supply one directory path.", call. = FALSE)
  }
  if (create) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.copy_engine <- function(destination) {
  source <- .engine_dir()
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(source, full.names = TRUE)
  for (file in files) {
    target <- file.path(destination, basename(file))
    if (file.exists(target)) {
      if (!identical(unname(tools::md5sum(file)), unname(tools::md5sum(target)))) {
        stop("Workspace engine differs: ", basename(file), call. = FALSE)
      }
    } else if (!file.copy(file, target)) {
      stop("Could not copy engine file: ", basename(file), call. = FALSE)
    }
  }
  invisible(destination)
}

.config_args <- function(config) {
  fields <- c("reps", "K_values", "folds", "base_seed", "mc_draws_ml",
              "gcomp_boot", "reference_draws", "canonical_estimators",
              "lambda_laws", "classifier_names")
  flags <- c("reps", "K-values", "folds", "base-seed", "mc-draws-ml",
             "gcomp-boot", "reference-draws", "canonical-estimators",
             "lambda-laws", "classifier-names")
  unlist(Map(function(field, flag) {
    c(paste0("--", flag), paste(config[[field]], collapse = ","))
  }, fields, flags), use.names = FALSE)
}

.dependencies <- function(packages) {
  present <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (!all(present)) {
    stop("Install required packages: ", paste(packages[!present], collapse = ", "),
         call. = FALSE)
  }
}

.execute <- function(engine, script, args = character(), log) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(engine)
  # Select the executable belonging to the active R installation.
  executable <- if (.Platform$OS.type == "windows") {
    file.path(R.home("bin"), "Rterm.exe")
  } else file.path(R.home("bin"), "exec", "R")
  if (!file.exists(executable)) {
    stop("R executable is unavailable: ", executable, call. = FALSE)
  }
  variables <- c(
    R_HOME = R.home(),
    R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep),
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1",
    R_DATATABLE_NUM_THREADS = "1"
  )
  old_env <- Sys.getenv(names(variables), unset = NA_character_)
  on.exit({
    Sys.unsetenv(names(old_env)[is.na(old_env)])
    restore <- old_env[!is.na(old_env)]
    if (length(restore)) do.call(Sys.setenv, as.list(restore))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(variables))
  dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  status <- system2(executable, shQuote(c(
    "--vanilla", "--slave", "-f", script, "--args", args
  )), stdout = log, stderr = log)
  if (status != 0L) {
    lines <- if (file.exists(log)) tail(readLines(log, warn = FALSE), 16L) else ""
    stop("Engine execution failed. Log: ", log, "\n",
         paste(lines, collapse = "\n"), call. = FALSE)
  }
  invisible(log)
}

.load_workspace <- function(path) {
  path <- .absolute_dir(path)
  config <- .validate_config(readRDS(file.path(path, "simulation_config.rds")))
  .copy_engine(file.path(path, "engine"))
  list(path = path, config = config)
}

