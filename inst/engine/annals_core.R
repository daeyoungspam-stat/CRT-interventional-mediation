# Shared mathematical and computational core for the Annals-level simulation.
# A cluster is the independent sampling unit throughout this project.

options(stringsAsFactors = FALSE)

ap_default_config <- function() {
  list(
    version = "annals-v8-sl-glm-ranger-dgp9",
    reps = 2000L,
    K_values = c(50L, 75L, 100L, 150L),
    folds = 5L,
    base_seed = 202608220L,
    ci_alpha = 0.05,
    ci_df_adjust = 3L,
    prob_lower = 0.025,
    prob_upper = 0.975,
    ratio_lower = 0.05,
    ratio_upper = 20,
    mc_draws_ml = 80L,
    gcomp_boot = 200L,
    reference_draws = 10L,
    canonical_estimators = c("par", "ml", "gcomp"),
    canonical_ml_library = c("SL.glm", "SL.ranger"),
    lowdim_s_library = c("SL.glm", "SL.ranger"),
    lambda_laws = c("gaussian", "resid_boot"),
    classifier_names = c("glmnet", "ranger", "sl"),
    dgp = list(
      n_min = 10L,
      n_max = 20L,
      treatment_prob = 0.5,
      v_n_coef = 0.20,
      x_v_coef = 0.55,
      x_n_coef = 0.15,
      x_sd = 1.00,
      z_intercept = -0.45,
      z_n = 0.35,
      z_v = 0.65,
      z_xbar = 0.55,
      z_a = 0.152,
      z_a_n = 0.072,
      z_a_v = 0.000,
      z_a_xbar = 0.000,
      m_intercept = -0.30,
      m_n = -0.30,
      m_v = 0.65,
      m_xbar = -0.55,
      m_a = 0.192,
      m_a_n = 0.000,
      m_a_v = 0.000,
      m_a_xbar = 0.000,
      sigma_z = 0.85,
      sigma_m_error = 0.90,
      gamma_zm = 0.10,
      process_rho = 0.06,
      eta_process_scale = 1.80,
      eta_zm_treatment_modifier = -0.15,
      sigma_y = 1.20,
      outcome_rho = 0.08
    )
  )
}

ap_stop <- function(...) stop(..., call. = FALSE)

ap_assert <- function(condition, message) {
  if (!isTRUE(condition)) ap_stop(message)
  invisible(TRUE)
}

ap_clip <- function(x, lower, upper) {
  pmin(pmax(as.numeric(x), lower), upper)
}

ap_clip_prob <- function(x, config) {
  ap_clip(x, config$prob_lower, config$prob_upper)
}

ap_clip_ratio <- function(x, config) {
  ap_clip(x, config$ratio_lower, config$ratio_upper)
}

ap_ml_integration_draws <- function(config, K) {
  # Let numerical integration accuracy increase with the independent sample
  # size so that Monte Carlo integration does not remain in the first-order
  # limit.  config$mc_draws_ml is the finite-sample lower bound.
  max(as.integer(config$mc_draws_ml), as.integer(K))
}

ap_hash_seed <- function(base_seed, ...) {
  key <- paste(..., collapse = "|")
  hash <- as.double(base_seed %% 2147483647)
  for (value in utf8ToInt(key)) hash <- (hash * 131 + value) %% 2147483647
  seed <- as.integer(hash)
  if (!is.finite(seed) || seed < 1L) seed <- 1L
  seed
}

ap_set_seed <- function(base_seed, ...) {
  seed <- ap_hash_seed(base_seed, ...)
  set.seed(seed)
  invisible(seed)
}

ap_parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (startsWith(token, "--")) {
      token <- sub("^--", "", token)
      if (grepl("=", token, fixed = TRUE)) {
        pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
        out[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
      } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        out[[token]] <- args[[i + 1L]]
        i <- i + 1L
      } else {
        out[[token]] <- "TRUE"
      }
    }
    i <- i + 1L
  }
  out
}

ap_arg <- function(args, name, default = NULL) {
  if (is.null(args[[name]])) default else args[[name]]
}

ap_config_signature <- function(config) {
  payload <- paste(capture.output(dput(config)), collapse = "|")
  hash <- 0
  for (value in utf8ToInt(payload)) hash <- (hash * 131 + value) %% 2147483647
  sprintf("%s-%010d", config$version, as.integer(hash))
}

ap_estimator_code_md5 <- function(root = ".") {
  filenames <- c("annals_core.R", "annals_canonical.R", "annals_lowdim.R")
  paths <- file.path(root, filenames)
  ap_assert(all(file.exists(paths)), paste(
    "Missing estimator source file(s):",
    paste(filenames[!file.exists(paths)], collapse = ", ")
  ))
  hashes <- unname(tools::md5sum(paths))
  stats::setNames(hashes, filenames)
}

ap_estimator_code_signature <- function(root = ".") {
  hashes <- ap_estimator_code_md5(root)
  payload <- paste(names(hashes), hashes, collapse = "|")
  hash <- 0
  for (value in utf8ToInt(payload)) hash <- (hash * 131 + value) %% 2147483647
  sprintf("estimator-code-%010d", as.integer(hash))
}

ap_expected_rows_per_rep <- function(study, config) {
  targets <- 6L
  if (study == "canonical") {
    return(length(ap_canonical_scenarios()) *
      length(config$canonical_estimators) * targets)
  }
  if (study == "lowdim") {
    return(length(ap_lowdim_scenarios()) * length(config$lambda_laws) *
      length(config$classifier_names) * targets)
  }
  ap_stop("Unknown study: ", study)
}

ap_result_path <- function(output_dir, study, K, rep_id) {
  file.path(
    output_dir, "parts", study, paste0("K", K),
    sprintf("rep_%04d.rds", rep_id)
  )
}

ap_n_center <- function(N) (as.numeric(N) - 15) / 5

ap_draw_exchangeable <- function(n, sigma, rho) {
  ap_assert(n >= 2L, "Cluster size must be at least two.")
  sigma * (sqrt(rho) * rnorm(1L) + sqrt(1 - rho) * rnorm(n))
}

ap_peer_mean <- function(x) {
  n <- length(x)
  (sum(x) - x) / (n - 1)
}

ap_process_mean <- function(data, arm, process = c("z", "m"), config) {
  process <- match.arg(process)
  d <- config$dgp
  A <- rep(as.numeric(arm), nrow(data))
  n <- ap_n_center(data$N)
  if (process == "z") {
    return(d$z_intercept + d$z_n * n + d$z_v * data$V +
      d$z_xbar * data$X_mean + A * (d$z_a + d$z_a_n * n +
      d$z_a_v * data$V + d$z_a_xbar * data$X_mean))
  }
  d$m_intercept + d$m_n * n + d$m_v * data$V +
    d$m_xbar * data$X_mean + A * (d$m_a + d$m_a_n * n +
    d$m_a_v * data$V + d$m_a_xbar * data$X_mean)
}

ap_eta_coefficients <- function(data, arm, config = NULL) {
  A <- rep(as.numeric(arm), nrow(data))
  scale <- if (is.null(config)) 1.80 else config$dgp$eta_process_scale
  zm_treatment_modifier <- if (is.null(config)) -0.15 else
    config$dgp$eta_zm_treatment_modifier
  list(
    base = 0.25 + 0.70 * A + 0.18 * ap_n_center(data$N) +
      0.22 * data$V + 0.15 * data$X + 0.10 * data$X_mean +
      0.12 * A * data$V - 0.10 * A * data$X_mean,
    z = scale * (0.50 + 0.12 * data$V + 0.12 * A),
    zp = scale * (0.30 + 0.10 * A),
    m = scale * (1.10 + 0.14 * data$V + 0.20 * A),
    mp = scale * (0.40 + 0.12 * A),
    z2 = rep(scale * 0.06, nrow(data)),
    zp2 = rep(0, nrow(data)),
    m2 = rep(scale * 0.10, nrow(data)),
    mp2 = rep(0, nrow(data)),
    z_zp = rep(0, nrow(data)),
    m_mp = rep(scale * 0.25, nrow(data)),
    z_m = scale * (0.70 + zm_treatment_modifier * A),
    zp_mp = rep(0, nrow(data)),
    z_mp = rep(0, nrow(data)),
    zp_m = rep(0, nrow(data))
  )
}

ap_eta_true <- function(data, arm = NULL, config = NULL) {
  if (is.null(arm)) arm <- data$A
  ap_assert(length(unique(arm)) == 1L || length(arm) == nrow(data),
            "arm must be scalar or row-specific.")
  if (length(unique(arm)) > 1L) {
    answer <- numeric(nrow(data))
    for (a in 0:1) {
      index <- which(arm == a)
      answer[index] <- ap_eta_true(data[index, , drop = FALSE], a, config)
    }
    return(answer)
  }
  b <- ap_eta_coefficients(data, arm[[1L]], config)
  with(data,
    b$base + b$z * Z + b$zp * Z_peer + b$m * M + b$mp * M_peer +
      b$z2 * Z^2 + b$zp2 * Z_peer^2 + b$m2 * M^2 +
      b$mp2 * M_peer^2 + b$z_zp * Z * Z_peer +
      b$m_mp * M * M_peer + b$z_m * Z * M +
      b$zp_mp * Z_peer * M_peer + b$z_mp * Z * M_peer +
      b$zp_m * Z_peer * M)
}

ap_generate_data <- function(K, config, rep_id, include_y = TRUE) {
  d <- config$dgp
  ap_set_seed(config$base_seed, "dgp", K, rep_id)
  clusters <- vector("list", K)
  for (id in seq_len(K)) {
    N <- sample(seq.int(d$n_min, d$n_max), 1L)
    nc <- ap_n_center(N)
    V <- rnorm(1L, mean = d$v_n_coef * nc, sd = 1)
    X <- d$x_v_coef * V + d$x_n_coef * nc + rnorm(N, sd = d$x_sd)
    X_mean <- mean(X)
    A <- rbinom(1L, 1L, d$treatment_prob)
    base <- data.frame(
      id = id, j = seq_len(N), N = N, A = A, V = V,
      X = X, X_mean = X_mean
    )
    mu_z <- ap_process_mean(base, A, "z", config)
    z_error <- ap_draw_exchangeable(N, d$sigma_z, d$process_rho)
    Z <- mu_z + z_error
    mu_m <- ap_process_mean(base, A, "m", config)
    m_error <- ap_draw_exchangeable(N, d$sigma_m_error, d$process_rho)
    M <- mu_m + d$gamma_zm * z_error + m_error
    cluster <- transform(base,
      Z = Z, Z_peer = ap_peer_mean(Z), M = M, M_peer = ap_peer_mean(M)
    )
    cluster$eta0 <- ap_eta_true(cluster, A, config)
    if (include_y) {
      cluster$Y <- cluster$eta0 +
        ap_draw_exchangeable(N, d$sigma_y, d$outcome_rho)
    } else {
      cluster$Y <- cluster$eta0
    }
    clusters[[id]] <- cluster
  }
  data <- do.call(rbind, clusters)
  rownames(data) <- NULL
  data
}

ap_cluster_frame <- function(data) {
  first <- !duplicated(data$id)
  data.frame(
    id = data$id[first], N = data$N[first], A = data$A[first],
    V = data$V[first], X_mean = data$X_mean[first]
  )
}

ap_make_folds <- function(data, V, seed) {
  cf <- ap_cluster_frame(data)
  fold <- integer(nrow(cf))
  set.seed(seed)
  for (a in 0:1) {
    index <- which(cf$A == a)
    index <- sample(index, length(index))
    fold[index] <- rep(seq_len(V), length.out = length(index))
  }
  setNames(fold, cf$id)
}

ap_eq_log_density <- function(y, mean, sigma, rho) {
  y <- as.numeric(y)
  mean <- as.numeric(mean)
  n <- length(y)
  ap_assert(length(mean) == n, "Mean and observation dimensions differ.")
  e <- (y - mean) / sigma
  denominator <- 1 + (n - 1) * rho
  logdet <- (n - 1) * log1p(-rho) + log(denominator)
  quad <- (sum(e^2) - rho * sum(e)^2 / denominator) / (1 - rho)
  -n * log(sigma) - n * log(2 * pi) / 2 - logdet / 2 - quad / 2
}

ap_eq_leave_log_density <- function(y, mean, sigma, rho) {
  y <- as.numeric(y)
  mean <- as.numeric(mean)
  n <- length(y)
  ap_assert(n >= 2L, "Leave-one density requires N >= 2.")
  e <- (y - mean) / sigma
  q <- n - 1L
  sums <- sum(e) - e
  sums2 <- sum(e^2) - e^2
  denominator <- 1 + (q - 1) * rho
  logdet <- (q - 1) * log1p(-rho) + log(denominator)
  quad <- (sums2 - rho * sums^2 / denominator) / (1 - rho)
  -q * log(sigma) - q * log(2 * pi) / 2 - logdet / 2 - quad / 2
}

ap_eq_marginal_conditional <- function(y, mean, sigma, rho) {
  n <- length(y)
  residual <- y - mean
  denominator <- 1 + (n - 2) * rho
  conditional_mean <- mean + rho * (sum(residual) - residual) / denominator
  variance_factor <- 1 - rho^2 * (n - 1) / denominator
  list(
    marginal_log = dnorm(y, mean, sigma, log = TRUE),
    conditional_log = dnorm(y, conditional_mean,
      sigma * sqrt(variance_factor), log = TRUE),
    conditional_mean = conditional_mean,
    conditional_sd = sigma * sqrt(variance_factor)
  )
}

ap_process_features <- function(data, arm = NULL, correct = TRUE,
                                conditional = FALSE) {
  A <- if (is.null(arm)) data$A else rep(as.numeric(arm), nrow(data))
  nc <- ap_n_center(data$N)
  if (isTRUE(correct)) {
    out <- data.frame(
      A = A, n = nc, V = data$V, X_mean = data$X_mean,
      A_n = A * nc
    )
    if (conditional) out <- cbind(data.frame(Z = data$Z), out)
    return(out)
  }
  # Reduced working laws retain the defining conditioning arguments and cluster
  # size, but omit baseline covariates and treatment-by-size modification.
  out <- data.frame(A = A, n = nc)
  if (conditional) out <- cbind(data.frame(Z = data$Z), out)
  out
}

ap_eta_features <- function(data, correct = TRUE) {
  if (!isTRUE(correct)) {
    return(data.frame(
      n = ap_n_center(data$N), Z = data$Z, Z_peer = data$Z_peer,
      M = data$M, M_peer = data$M_peer,
      Z2 = data$Z^2, M2 = data$M^2,
      M_M_peer = data$M * data$M_peer,
      Z_M = data$Z * data$M
    ))
  }
  data.frame(
    n = ap_n_center(data$N), V = data$V, X = data$X,
    X_mean = data$X_mean, Z = data$Z, Z_peer = data$Z_peer,
    M = data$M, M_peer = data$M_peer,
    V_Z = data$V * data$Z,
    V_M = data$V * data$M,
    Z2 = data$Z^2,
    M2 = data$M^2,
    M_M_peer = data$M * data$M_peer,
    Z_M = data$Z * data$M
  )
}

ap_model_frame <- function(X) {
  X <- as.data.frame(X)
  X[] <- lapply(X, function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- 0
    x
  })
  X
}

ap_full_rank_model_frame <- function(X, tolerance = 1e-10) {
  X <- ap_model_frame(X)
  design <- cbind(`(Intercept)` = 1, as.matrix(X))
  decomposition <- qr(design, tol = tolerance, LAPACK = FALSE)
  retained <- decomposition$pivot[seq_len(decomposition$rank)]
  retained <- retained[retained != 1L] - 1L
  X[, retained, drop = FALSE]
}

ap_fit_base_model <- function(Y, X, family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  X <- ap_model_frame(X)
  design <- cbind(`(Intercept)` = 1, as.matrix(X))
  if (family == "gaussian") {
    fit <- lm.fit(design, as.numeric(Y))
    beta <- fit$coefficients
    beta[!is.finite(beta)] <- 0
    return(list(type = "base_lm", coefficients = beta,
                names = colnames(design), family = family))
  }
  fit <- suppressWarnings(glm.fit(design, as.numeric(Y), family = binomial()))
  beta <- fit$coefficients
  beta[!is.finite(beta)] <- 0
  list(type = "base_glm", coefficients = beta,
       names = colnames(design), family = family)
}

ap_predict_base_model <- function(model, X, linear = FALSE) {
  X <- ap_model_frame(X)
  design <- cbind(`(Intercept)` = 1, as.matrix(X))
  missing <- setdiff(model$names, colnames(design))
  if (length(missing)) {
    design <- cbind(design, matrix(0, nrow(design), length(missing),
      dimnames = list(NULL, missing)))
  }
  design <- design[, model$names, drop = FALSE]
  eta <- drop(design %*% model$coefficients)
  if (model$family == "binomial" && !linear) return(plogis(eta))
  eta
}

predict.ap_SL_glm_vridge <- function(object, newdata, ...) {
  newdata <- ap_model_frame(newdata)
  missing <- setdiff(object$names, names(newdata))
  for (name in missing) newdata[[name]] <- 0
  newdata <- newdata[, object$names, drop = FALSE]
  as.numeric(predict(
    object$object, newx = as.matrix(newdata), s = object$lambda,
    type = "response"
  ))
}

ap_make_SL_glm_vridge <- function(lambda, audit) {
  force(lambda)
  force(audit)
  function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
    X <- ap_model_frame(X)
    newX <- ap_model_frame(newX)
    if (!identical(family$family, "binomial")) {
      return(SuperLearner::SL.glm(
        Y, X, newX, family, obsWeights, model = model, ...
      ))
    }

    # Exact linear dependencies do not change the logistic-regression span.
    # Remove them inside the GLM wrapper only, leaving the ranger learner's
    # prespecified feature set untouched.
    glm_X <- ap_full_rank_model_frame(X)
    glm_newX <- newX[, names(glm_X), drop = FALSE]
    warning_messages <- character(0)
    ordinary <- tryCatch(
      withCallingHandlers({
        fit <- stats::glm(
          Y ~ ., data = glm_X, family = family, weights = obsWeights,
          model = model
        )
        prediction <- stats::predict(
          fit, newdata = glm_newX, type = "response"
        )
        list(fit = fit, prediction = as.numeric(prediction))
      }, warning = function(warning) {
        warning_messages <<- c(warning_messages, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }),
      error = function(error) error
    )
    trigger_pattern <- paste(
      c("did not converge", "numerically 0 or 1", "stopped at boundary"),
      collapse = "|"
    )
    boundary_tolerance <- 1e-8
    use_ridge <- inherits(ordinary, "error") ||
      length(grep(trigger_pattern, warning_messages, ignore.case = TRUE)) > 0L ||
      !isTRUE(ordinary$fit$converged) ||
      any(!is.finite(ordinary$fit$coefficients)) ||
      any(!is.finite(ordinary$prediction)) ||
      any(ordinary$prediction <= boundary_tolerance |
          ordinary$prediction >= 1 - boundary_tolerance)

    if (!use_ridge) {
      fit <- list(object = ordinary$fit)
      class(fit) <- "SL.glm"
      return(list(pred = ordinary$prediction, fit = fit))
    }

    audit$count <- audit$count + 1L
    ridge <- glmnet::glmnet(
      x = as.matrix(glm_X), y = as.numeric(Y), family = "binomial",
      weights = obsWeights, alpha = 0, lambda = lambda,
      standardize = TRUE, intercept = TRUE,
      control = list(maxit = 1000000L)
    )
    prediction <- as.numeric(predict(
      ridge, newx = as.matrix(glm_newX), s = lambda, type = "response"
    ))
    ap_assert(all(is.finite(prediction)),
              "Vanishing-ridge fallback produced non-finite predictions.")
    fit <- list(object = ridge, lambda = lambda, names = names(glm_X))
    class(fit) <- "ap_SL_glm_vridge"
    list(pred = prediction, fit = fit)
  }
}

ap_fit_model <- function(Y, X, family = "gaussian", learner = "par",
                         id = NULL, seed = 1L, library = NULL) {
  if (learner == "par") return(ap_fit_base_model(Y, X, family))
  ap_assert(requireNamespace("SuperLearner", quietly = TRUE),
            "SuperLearner is required for learner != 'par'.")
  mapping <- list(
    ml = library,
    glmnet = "SL.glmnet",
    ranger = "SL.ranger",
    sl = c("SL.glm", "SL.ranger")
  )
  sl_library <- mapping[[learner]]
  ap_assert(length(sl_library) > 0L, paste("Unknown learner:", learner))
  available <- vapply(sl_library, function(name) {
    switch(name,
      SL.glm = TRUE,
      SL.glmnet = requireNamespace("glmnet", quietly = TRUE),
      SL.ranger = requireNamespace("ranger", quietly = TRUE),
      FALSE)
  }, logical(1))
  sl_library <- sl_library[available]
  ap_assert(length(sl_library) > 0L, "No requested SuperLearner wrapper is available.")
  X <- ap_model_frame(X)
  set.seed(seed)
  unique_units <- if (is.null(id)) length(Y) else length(unique(id))
  ridge_audit <- new.env(parent = emptyenv())
  ridge_audit$count <- 0L
  sl_environment <- asNamespace("SuperLearner")
  ridge_lambda <- NA_real_
  if ("SL.glm" %in% sl_library) {
    ap_assert(requireNamespace("glmnet", quietly = TRUE),
              "glmnet is required for the prespecified SL.glm ridge fallback.")
    ridge_lambda <- 1 / unique_units
    sl_environment <- new.env(parent = asNamespace("SuperLearner"))
    sl_environment$SL.glm <- ap_make_SL_glm_vridge(ridge_lambda, ridge_audit)
  }
  warning_messages <- character(0)
  fit <- withCallingHandlers(
    SuperLearner::SuperLearner(
      Y = as.numeric(Y), X = X,
      family = if (family == "gaussian") gaussian() else binomial(),
      id = id, SL.library = sl_library,
      cvControl = list(V = max(2L, min(5L, unique_units)), shuffle = TRUE),
      control = list(saveFitLibrary = TRUE, trimLogit = 0.01),
      env = sl_environment
    ),
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  degenerate <- !length(fit$coef) || any(!is.finite(fit$coef)) ||
    sum(abs(fit$coef)) <= 1e-12
  if (degenerate) {
    fallback <- ap_fit_base_model(Y, X, family)
    fallback$fallback <- "degenerate_superlearner"
    fallback$superlearner_warnings <- unique(warning_messages)
    return(fallback)
  }
  if (length(warning_messages)) {
    warning(paste(unique(warning_messages), collapse = "; "), call. = FALSE)
  }
  list(type = "SuperLearner", fit = fit, names = names(X),
       family = family, library = sl_library,
       ridge_fallback_count = ridge_audit$count,
       ridge_lambda = ridge_lambda)
}

ap_predict_model <- function(model, X, linear = FALSE) {
  if (model$type %in% c("base_lm", "base_glm")) {
    return(ap_predict_base_model(model, X, linear))
  }
  X <- ap_model_frame(X)
  missing <- setdiff(model$names, names(X))
  for (name in missing) X[[name]] <- 0
  X <- X[, model$names, drop = FALSE]
  as.numeric(predict(model$fit, newdata = X, onlySL = TRUE)$pred)
}

ap_estimate_eq_covariance <- function(residual, id) {
  residual <- as.numeric(residual)
  sigma2 <- mean(residual^2)
  sigma2 <- max(sigma2, 1e-8)
  split_residual <- split(residual, id)
  numerator <- sum(vapply(split_residual, function(x) sum(x)^2 - sum(x^2), 0.0))
  denominator <- sigma2 * sum(vapply(split_residual,
    function(x) length(x) * (length(x) - 1), 0.0))
  rho <- if (denominator > 0) numerator / denominator else 0
  max_n <- max(vapply(split_residual, length, 1L))
  lower <- -0.95 / (max_n - 1)
  rho <- ap_clip(rho, max(lower, -0.05), 0.80)
  list(sigma = sqrt(sigma2), rho = rho)
}

ap_fit_process <- function(data, outcome = c("Z", "M"), correct = TRUE,
                           conditional = FALSE, learner = "par", seed = 1L,
                           library = NULL) {
  outcome <- match.arg(outcome)
  X <- ap_process_features(data, correct = correct, conditional = conditional)
  model <- ap_fit_model(data[[outcome]], X, "gaussian", learner,
    data$id, seed, library)
  mean <- ap_predict_model(model, X)
  covariance <- ap_estimate_eq_covariance(data[[outcome]] - mean, data$id)
  list(
    outcome = outcome, correct = correct, conditional = conditional,
    learner = learner, model = model, sigma = covariance$sigma,
    rho = covariance$rho
  )
}

ap_predict_process <- function(fit, data, arm) {
  X <- ap_process_features(data, arm, fit$correct, fit$conditional)
  ap_predict_model(fit$model, X)
}

ap_fit_eta <- function(data, arm, correct = TRUE, learner = "par",
                       seed = 1L, library = NULL) {
  # Both working laws are arm-specific.  The reduced law retains the own and
  # peer post-treatment arguments and their nonlinear terms but omits baseline
  # covariates and baseline effect modification; causal contrasts therefore
  # remain nondegenerate under deliberate misspecification.
  subset <- data[data$A == arm, , drop = FALSE]
  model <- ap_fit_model(subset$Y, ap_eta_features(subset, correct),
    "gaussian", learner, subset$id, seed, library)
  list(arm = arm, correct = correct, learner = learner, model = model)
}

ap_predict_eta <- function(fit, data) {
  if (identical(fit$type, "truth_eta")) {
    return(ap_eta_true(data, fit$arm, fit$config))
  }
  ap_predict_model(fit$model, ap_eta_features(data, fit$correct))
}

ap_true_eta_object <- function(arm, config) {
  list(type = "truth_eta", arm = as.numeric(arm), correct = TRUE,
       learner = "par", config = config)
}

ap_empty_moments <- function(n) {
  zero <- numeric(n)
  list(
    z = zero, zp = zero, m = zero, mp = zero,
    var_z = zero, var_zp = zero, var_m = zero, var_mp = zero,
    cov_z_zp = zero, cov_m_mp = zero,
    cov_z_m = zero, cov_z_mp = zero,
    cov_zp_m = zero, cov_zp_mp = zero
  )
}

ap_eq_pair_moments <- function(mean, sigma, rho, N) {
  q <- N - 1
  list(
    own_mean = mean,
    peer_mean = mean,
    own_var = rep(sigma^2, length(N)),
    peer_var = sigma^2 * (rho + (1 - rho) / q),
    own_peer_cov = rep(sigma^2 * rho, length(N))
  )
}

ap_true_process_object <- function(process = c("phi", "kappa", "xi"), config) {
  process <- match.arg(process)
  d <- config$dgp
  if (process == "phi") {
    return(list(type = "truth", outcome = "Z", conditional = FALSE,
      sigma = d$sigma_z, rho = d$process_rho, config = config))
  }
  if (process == "kappa") {
    sigma2 <- d$gamma_zm^2 * d$sigma_z^2 + d$sigma_m_error^2
    return(list(type = "truth", outcome = "M", conditional = FALSE,
      sigma = sqrt(sigma2), rho = d$process_rho, config = config))
  }
  list(type = "truth", outcome = "M", conditional = TRUE,
    sigma = d$sigma_m_error, rho = d$process_rho, config = config)
}

ap_process_prediction <- function(process, data, arm) {
  if (identical(process$type, "derived_kappa")) {
    at_marginal_mean <- data
    at_marginal_mean$Z <- ap_process_prediction(process$phi, data, arm)
    return(ap_process_prediction(process$xi, at_marginal_mean, arm))
  }
  if (!identical(process$type, "truth")) return(ap_predict_process(process, data, arm))
  if (process$outcome == "Z") return(ap_process_mean(data, arm, "z", process$config))
  mu_m <- ap_process_mean(data, arm, "m", process$config)
  if (!isTRUE(process$conditional)) return(mu_m)
  mu_z <- ap_process_mean(data, arm, "z", process$config)
  mu_m + process$config$dgp$gamma_zm * (data$Z - mu_z)
}

ap_derive_compatible_kappa <- function(phi, xi, data) {
  # Under the fitted linear-Gaussian joint law, integrate the conditional
  # M|Z law against the fitted marginal Z law.  This makes the fitted kappa
  # exactly compatible with the fitted phi and xi, including its covariance.
  probe <- data
  baseline <- ap_predict_process(xi, probe, 0)
  probe$Z <- probe$Z + 1
  increment <- ap_predict_process(xi, probe, 0) - baseline
  gamma <- mean(increment)
  ap_assert(max(abs(increment - gamma)) < 1e-8,
            "The correct conditional mediator model is not linear in Z.")
  sigma2 <- gamma^2 * phi$sigma^2 + xi$sigma^2
  rho <- (gamma^2 * phi$sigma^2 * phi$rho +
    xi$sigma^2 * xi$rho) / sigma2
  list(
    type = "derived_kappa", outcome = "M", conditional = FALSE,
    correct = TRUE, learner = xi$learner, phi = phi, xi = xi,
    gamma = gamma, sigma = sqrt(sigma2), rho = rho
  )
}

ap_eta_feature_expectation <- function(data, moments, correct = TRUE) {
  if (!isTRUE(correct)) {
    return(data.frame(
      n = ap_n_center(data$N), Z = moments$z, Z_peer = moments$zp,
      M = moments$m, M_peer = moments$mp,
      Z2 = moments$z^2 + moments$var_z,
      M2 = moments$m^2 + moments$var_m,
      M_M_peer = moments$m * moments$mp + moments$cov_m_mp,
      Z_M = moments$z * moments$m + moments$cov_z_m
    ))
  }
  data.frame(
    n = ap_n_center(data$N), V = data$V, X = data$X,
    X_mean = data$X_mean, Z = moments$z, Z_peer = moments$zp,
    M = moments$m, M_peer = moments$mp,
    V_Z = data$V * moments$z,
    V_M = data$V * moments$m,
    Z2 = moments$z^2 + moments$var_z,
    M2 = moments$m^2 + moments$var_m,
    M_M_peer = moments$m * moments$mp + moments$cov_m_mp,
    Z_M = moments$z * moments$m + moments$cov_z_m
  )
}

ap_predict_eta_moment <- function(fit, data, moments) {
  if (identical(fit$type, "truth_eta")) {
    return(ap_true_eta_moment(data, fit$arm, moments, fit$config))
  }
  ap_assert(fit$learner == "par",
            "Exact moment integration is available only for parametric eta fits.")
  X <- ap_eta_feature_expectation(data, moments, fit$correct)
  ap_predict_model(fit$model, X)
}

ap_true_eta_moment <- function(data, arm, moments, config = NULL) {
  b <- ap_eta_coefficients(data, arm, config)
  b$base + b$z * moments$z + b$zp * moments$zp +
    b$m * moments$m + b$mp * moments$mp +
    b$z2 * (moments$z^2 + moments$var_z) +
    b$zp2 * (moments$zp^2 + moments$var_zp) +
    b$m2 * (moments$m^2 + moments$var_m) +
    b$mp2 * (moments$mp^2 + moments$var_mp) +
    b$z_zp * (moments$z * moments$zp + moments$cov_z_zp) +
    b$m_mp * (moments$m * moments$mp + moments$cov_m_mp) +
    b$z_m * (moments$z * moments$m + moments$cov_z_m) +
    b$zp_mp * (moments$zp * moments$mp + moments$cov_zp_mp) +
    b$z_mp * (moments$z * moments$mp + moments$cov_z_mp) +
    b$zp_m * (moments$zp * moments$m + moments$cov_zp_m)
}

ap_process_pair <- function(process, data, arm) {
  mean <- ap_process_prediction(process, data, arm)
  ap_eq_pair_moments(mean, process$sigma, process$rho, data$N)
}

ap_moments_theta <- function(data, phi, kappa, a, a_prime,
                             keep_z = FALSE, keep_m = FALSE) {
  out <- ap_empty_moments(nrow(data))
  if (keep_z) {
    out$z <- data$Z
    out$zp <- data$Z_peer
  } else {
    z <- ap_process_pair(phi, data, a)
    out$z <- z$own_mean
    out$zp <- z$peer_mean
    out$var_z <- z$own_var
    out$var_zp <- z$peer_var
    out$cov_z_zp <- z$own_peer_cov
  }
  if (keep_m) {
    out$m <- data$M
    out$mp <- data$M_peer
  } else {
    m <- ap_process_pair(kappa, data, a_prime)
    out$m <- m$own_mean
    out$mp <- m$peer_mean
    out$var_m <- m$own_var
    out$var_mp <- m$peer_var
    out$cov_m_mp <- m$own_peer_cov
  }
  out
}

ap_moments_tau <- function(data, phi, kappa, keep_z = FALSE,
                           keep_own_m = FALSE, keep_peer_m = FALSE) {
  out <- ap_empty_moments(nrow(data))
  if (keep_z) {
    out$z <- data$Z
    out$zp <- data$Z_peer
  } else {
    z <- ap_process_pair(phi, data, 1)
    out$z <- z$own_mean
    out$zp <- z$peer_mean
    out$var_z <- z$own_var
    out$var_zp <- z$peer_var
    out$cov_z_zp <- z$own_peer_cov
  }
  own <- ap_process_pair(kappa, data, 1)
  peer <- ap_process_pair(kappa, data, 0)
  if (keep_own_m) {
    out$m <- data$M
  } else {
    out$m <- own$own_mean
    out$var_m <- own$own_var
  }
  if (keep_peer_m) {
    out$mp <- data$M_peer
  } else {
    out$mp <- peer$peer_mean
    out$var_mp <- peer$peer_var
  }
  # The intervention defining tau draws own and peer mediators independently.
  out$cov_m_mp <- numeric(nrow(data))
  out
}

ap_partial_theta_par <- function(eta_fit, data, phi, kappa, a, a_prime) {
  list(
    u1 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_theta(data, phi, kappa, a, a_prime, keep_z = TRUE)),
    u2 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_theta(data, phi, kappa, a, a_prime, keep_m = TRUE)),
    u3 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_theta(data, phi, kappa, a, a_prime))
  )
}

ap_partial_tau_par <- function(eta_fit, data, phi, kappa) {
  list(
    r1 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_tau(data, phi, kappa, keep_z = TRUE)),
    r2 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_tau(data, phi, kappa, keep_own_m = TRUE)),
    r3 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_tau(data, phi, kappa, keep_peer_m = TRUE)),
    r4 = ap_predict_eta_moment(eta_fit, data,
      ap_moments_tau(data, phi, kappa))
  )
}

ap_gauss_hermite <- function(order = 20L) {
  order <- as.integer(order)
  J <- matrix(0, order, order)
  off <- sqrt(seq_len(order - 1L) / 2)
  J[cbind(seq_len(order - 1L), 2:order)] <- off
  J[cbind(2:order, seq_len(order - 1L))] <- off
  eig <- eigen(J, symmetric = TRUE)
  index <- order(eig$values)
  list(
    node = sqrt(2) * eig$values[index],
    weight = eig$vectors[1L, index]^2
  )
}

.ap_truth_cache <- new.env(parent = emptyenv())

ap_truth <- function(config, quadrature_order = 15L) {
  cache_key <- paste(ap_config_signature(config), as.integer(quadrature_order), sep = "|")
  if (exists(cache_key, envir = .ap_truth_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .ap_truth_cache, inherits = FALSE))
  }
  d <- config$dgp
  gh <- ap_gauss_hermite(quadrature_order)
  components_C <- setNames(numeric(4L), c("theta11", "theta10", "theta00", "tau"))
  components_I <- components_C
  mean_N <- mean(seq.int(d$n_min, d$n_max))
  for (N in seq.int(d$n_min, d$n_max)) {
    nc <- ap_n_center(N)
    for (iv in seq_along(gh$node)) {
      V <- d$v_n_coef * nc + gh$node[[iv]]
      for (ix in seq_along(gh$node)) {
        X_mean <- d$x_v_coef * V + d$x_n_coef * nc +
          d$x_sd * gh$node[[ix]] / sqrt(N)
        weight <- gh$weight[[iv]] * gh$weight[[ix]] /
          (d$n_max - d$n_min + 1)
        row <- data.frame(
          id = 1, j = 1, N = N, A = 0, V = V,
          X = X_mean, X_mean = X_mean,
          Z = 0, Z_peer = 0, M = 0, M_peer = 0
        )
        phi <- ap_true_process_object("phi", config)
        kappa <- ap_true_process_object("kappa", config)
        values <- c(
          theta11 = ap_true_eta_moment(row, 1,
            ap_moments_theta(row, phi, kappa, 1, 1), config),
          theta10 = ap_true_eta_moment(row, 1,
            ap_moments_theta(row, phi, kappa, 1, 0), config),
          theta00 = ap_true_eta_moment(row, 0,
            ap_moments_theta(row, phi, kappa, 0, 0), config),
          tau = ap_true_eta_moment(row, 1,
            ap_moments_tau(row, phi, kappa), config)
        )
        components_C <- components_C + weight * values
        components_I <- components_I + weight * N * values / mean_N
      }
    }
  }
  effects <- function(x) c(
    IDE = x[["theta10"]] - x[["theta00"]],
    IIME = x[["tau"]] - x[["theta10"]],
    ISME = x[["theta11"]] - x[["tau"]]
  )
  eC <- effects(components_C)
  eI <- effects(components_I)
  answer <- c(
    IDE_C = eC[["IDE"]], IIME_C = eC[["IIME"]], ISME_C = eC[["ISME"]],
    IDE_I = eI[["IDE"]], IIME_I = eI[["IIME"]], ISME_I = eI[["ISME"]]
  )
  assign(cache_key, answer, envir = .ap_truth_cache)
  answer
}

ap_canonical_scenarios <- function() {
  grid <- expand.grid(
    eta = c(TRUE, FALSE), phi = c(TRUE, FALSE),
    kappa = c(TRUE, FALSE), xi = c(TRUE, FALSE),
    KEEP.OUT.ATTRS = FALSE
  )
  grid <- grid[order(-grid$eta, -grid$phi, -grid$kappa, -grid$xi), ]
  output <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    x <- grid[i, ]
    id <- sprintf("C-%d%d%d%d", x$eta, x$phi, x$kappa, x$xi)
    theta_point <- (x$eta && x$phi) || (x$eta && x$kappa) ||
      (x$kappa && x$xi)
    tau_point <- x$kappa && (x$eta || x$xi)
    all_correct <- all(unlist(x))
    output[[i]] <- list(
      id = id, eta = x$eta, phi = x$phi, kappa = x$kappa, xi = x$xi,
      point_consistency_expected = c(
        IDE = theta_point,
        IIME = theta_point && tau_point,
        ISME = theta_point && tau_point
      ),
      asymptotic_linearity_expected = c(
        IDE = all_correct, IIME = all_correct, ISME = all_correct
      ),
      ci_validity_expected = c(
        IDE = all_correct, IIME = all_correct, ISME = all_correct
      ),
      gcomp_point_expected = c(
        IDE = x$eta && x$phi && x$kappa,
        IIME = x$eta && x$phi && x$kappa,
        ISME = x$eta && x$phi && x$kappa
      )
    )
  }
  names(output) <- vapply(output, `[[`, "", "id")
  output
}

ap_lowdim_scenarios <- function() {
  all_theta <- list(eta = TRUE, d = TRUE, s = TRUE)
  all_tau <- list(eta = TRUE, e = TRUE, s = TRUE)
  theta <- list(
    `LD-T2` = list(eta = TRUE, d = FALSE, s = FALSE),
    `LD-T3` = list(eta = FALSE, d = TRUE, s = TRUE),
    `LD-T4` = list(eta = FALSE, d = TRUE, s = FALSE),
    `LD-T5` = list(eta = FALSE, d = FALSE, s = TRUE),
    `LD-T6` = list(eta = FALSE, d = FALSE, s = FALSE)
  )
  tau <- list(
    `LD-U2` = list(eta = TRUE, e = FALSE, s = FALSE),
    `LD-U3` = list(eta = FALSE, e = TRUE, s = TRUE),
    `LD-U4` = list(eta = FALSE, e = TRUE, s = FALSE),
    `LD-U5` = list(eta = FALSE, e = FALSE, s = TRUE),
    `LD-U6` = list(eta = FALSE, e = FALSE, s = FALSE)
  )
  make <- function(id, theta_cfg, tau_cfg) {
    theta_cross_ok <- theta_cfg$eta || (theta_cfg$d && theta_cfg$s)
    # For theta(a,a), the treatment-classifier odds equal one identically.
    theta_same_ok <- theta_cfg$eta || theta_cfg$d
    tau_ok <- tau_cfg$eta || (tau_cfg$e && tau_cfg$s)
    all_correct <- id == "LD-AC"
    list(
      id = id, theta = theta_cfg, tau = tau_cfg,
      point_consistency_expected = c(
        IDE = theta_cross_ok && theta_same_ok,
        IIME = theta_cross_ok && tau_ok,
        ISME = theta_same_ok && tau_ok
      ),
      asymptotic_linearity_expected = c(
        IDE = all_correct, IIME = all_correct, ISME = all_correct
      ),
      ci_validity_expected = c(
        IDE = all_correct, IIME = all_correct, ISME = all_correct
      )
    )
  }
  out <- list(`LD-AC` = make("LD-AC", all_theta, all_tau))
  for (id in names(theta)) out[[id]] <- make(id, theta[[id]], all_tau)
  for (id in names(tau)) out[[id]] <- make(id, all_theta, tau[[id]])
  out
}

ap_effect_contrasts <- function() {
  list(
    IDE = c(theta10 = 1, theta00 = -1),
    IIME = c(tau = 1, theta10 = -1),
    ISME = c(theta11 = 1, tau = -1)
  )
}

ap_effect_rows <- function(components, N_by_id, config) {
  K <- length(N_by_id)
  truth <- ap_truth(config)
  rows <- list()
  index <- 1L
  for (effect in names(ap_effect_contrasts())) {
    contrast <- ap_effect_contrasts()[[effect]]
    raw_C <- Reduce(`+`, Map(function(name, value) {
      value * components[[name]]$C
    }, names(contrast), as.numeric(contrast)))
    raw_I <- Reduce(`+`, Map(function(name, value) {
      value * components[[name]]$I
    }, names(contrast), as.numeric(contrast)))
    point_C <- mean(raw_C)
    point_I <- mean(raw_I) / mean(N_by_id)
    if_C <- raw_C - point_C
    if_I <- (raw_I - N_by_id * point_I) / mean(N_by_id)
    se_C <- sqrt(mean(if_C^2) / K)
    se_I <- sqrt(mean(if_I^2) / K)
    df <- max(1L, K - config$ci_df_adjust)
    critical <- sqrt(K / df) * qt(1 - config$ci_alpha / 2, df)
    for (scale in c("C", "I")) {
      point <- if (scale == "C") point_C else point_I
      se <- if (scale == "C") se_C else se_I
      target <- paste0(effect, "_", scale)
      rows[[index]] <- data.frame(
        estimand = target, effect = effect,
        scale = if (scale == "C") "cluster" else "individual",
        point = point, se = se,
        ci_low = point - critical * se,
        ci_up = point + critical * se,
        true = unname(truth[[target]])
      )
      index <- index + 1L
    }
  }
  do.call(rbind, rows)
}

ap_finalize_rows <- function(rows, study, K, rep_id, scenario, method,
                             theory, diagnostics = list(),
                             ci_method = "cluster_if_t") {
  rows$study <- study
  rows$K <- K
  rows$rep <- rep_id
  rows$scenario <- scenario
  rows$method <- method
  rows$point_consistency_expected <-
    as.logical(theory$point_consistency_expected[rows$effect])
  rows$asymptotic_linearity_expected <-
    as.logical(theory$asymptotic_linearity_expected[rows$effect])
  rows$ci_validity_expected <-
    as.logical(theory$ci_validity_expected[rows$effect])
  rows$ci_method <- ci_method
  rows$cover <- as.numeric(rows$ci_low <= rows$true & rows$ci_up >= rows$true)
  rows$error <- rows$point - rows$true
  rows$ci_width <- rows$ci_up - rows$ci_low
  for (name in names(diagnostics)) rows[[name]] <- diagnostics[[name]]
  first <- c(
    "study", "K", "rep", "scenario", "method", "estimand", "effect",
    "scale", "point_consistency_expected", "asymptotic_linearity_expected",
    "ci_validity_expected", "ci_method"
  )
  rows[, c(first, setdiff(names(rows), first)), drop = FALSE]
}

ap_validate_config <- function(config) {
  ap_assert(identical(c(config$prob_lower, config$prob_upper), c(0.025, 0.975)),
            "Classifier truncation must remain exactly [0.025, 0.975].")
  ap_assert(identical(c(config$ratio_lower, config$ratio_upper), c(0.05, 20)),
            "Ratio truncation must remain exactly [0.05, 20].")
  ap_assert(config$dgp$treatment_prob == 0.5,
            "Current score implementation assumes randomized probability 0.5.")
  ap_assert(config$mc_draws_ml >= 1L,
            "The ML integration-draw lower bound must be positive.")
  ap_assert(identical(config$canonical_ml_library, c("SL.glm", "SL.ranger")),
            "Canonical ML must use exactly SL.glm and SL.ranger.")
  ap_assert(identical(config$lowdim_s_library, c("SL.glm", "SL.ranger")),
            "Low-dimensional Super Learner fits must use exactly SL.glm and SL.ranger.")
  ap_assert(config$reference_draws >= 1L,
            "The augmented-classification reference-draw count must be positive.")
  ap_assert(config$dgp$process_rho >= 0 && config$dgp$process_rho < 1,
            "Process rho must be in [0,1).")
  invisible(TRUE)
}
