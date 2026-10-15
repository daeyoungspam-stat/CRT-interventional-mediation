# Low-dimensional EIF based on augmented classification.

ap_lambda_features <- function(data) {
  data.frame(
    N = as.numeric(data$N), V = data$V,
    X_mean = data$X_mean
  )
}

ap_fit_lambda_reference <- function(train_arm, seed) {
  model <- ap_fit_model(
    train_arm$M, ap_lambda_features(train_arm), "gaussian", "par",
    train_arm$id, seed
  )
  mean <- ap_predict_model(model, ap_lambda_features(train_arm))
  residual <- train_arm$M - mean
  residual <- residual - mean(residual)
  list(
    model = model, residual = residual,
    sigma = max(sqrt(mean(residual^2)), 1e-6)
  )
}

ap_draw_lambda <- function(reference, cluster, law, seed) {
  mean <- ap_predict_model(reference$model, ap_lambda_features(cluster))
  set.seed(seed)
  if (law == "gaussian") {
    return(mean + rnorm(nrow(cluster), sd = reference$sigma))
  }
  if (law == "resid_boot") {
    return(mean + sample(reference$residual, nrow(cluster), replace = TRUE))
  }
  ap_stop("Unknown reference law: ", law)
}

ap_baseline_polynomial <- function(n, V, X_mean, multiplier = 1) {
  data.frame(
    base_0 = multiplier,
    base_n = multiplier * n,
    base_V = multiplier * V,
    base_X = multiplier * X_mean,
    base_n2 = multiplier * n^2,
    base_V2 = multiplier * V^2,
    base_X2 = multiplier * X_mean^2,
    base_nV = multiplier * n * V,
    base_nX = multiplier * n * X_mean,
    base_VX = multiplier * V * X_mean
  )
}

ap_d_full_features <- function(cluster, mediator = cluster$M, rho_basis) {
  z <- cluster$Z
  m <- as.numeric(mediator)
  N <- length(m)
  n <- ap_n_center(N)
  V <- cluster$V[[1L]]
  X_mean <- cluster$X_mean[[1L]]
  denominator <- 1 + (N - 1) * rho_basis
  sum_z <- sum(z)
  sum_m <- sum(m)
  base <- data.frame(
    N = N, n = n, V = V, X_mean = X_mean, log_D = log(denominator),
    sum_Z = sum_z, sum_M = sum_m,
    sum_Z2 = sum(z^2), sum_M2 = sum(m^2), sum_ZM = sum(z * m),
    sum_Z_sq_over_D = sum_z^2 / denominator,
    sum_M_sq_over_D = sum_m^2 / denominator,
    sum_ZM_over_D = sum_z * sum_m / denominator,
    sum_Z_over_D = sum_z / denominator,
    sum_M_over_D = sum_m / denominator,
    n_sum_Z_over_D = n * sum_z / denominator,
    V_sum_Z_over_D = V * sum_z / denominator,
    X_sum_Z_over_D = X_mean * sum_z / denominator,
    n_sum_M_over_D = n * sum_m / denominator,
    V_sum_M_over_D = V * sum_m / denominator,
    X_sum_M_over_D = X_mean * sum_m / denominator,
    n_sum_M = n * sum_m,
    V_sum_M = V * sum_m,
    X_sum_M = X_mean * sum_m
  )
  scaled <- ap_baseline_polynomial(n, V, X_mean, N / denominator)
  raw <- ap_baseline_polynomial(n, V, X_mean, N)
  names(raw) <- paste0("raw_", names(raw))
  cbind(base, scaled, raw)
}

ap_d_m_features <- function(cluster, mediator = cluster$M, rho_basis) {
  m <- as.numeric(mediator)
  N <- length(m)
  n <- ap_n_center(N)
  V <- cluster$V[[1L]]
  X_mean <- cluster$X_mean[[1L]]
  denominator <- 1 + (N - 1) * rho_basis
  sum_m <- sum(m)
  base <- data.frame(
    N = N, n = n, V = V, X_mean = X_mean, log_D = log(denominator),
    sum_M = sum_m, sum_M2 = sum(m^2),
    sum_M_sq_over_D = sum_m^2 / denominator,
    sum_M_over_D = sum_m / denominator,
    n_sum_M_over_D = n * sum_m / denominator,
    V_sum_M_over_D = V * sum_m / denominator,
    X_sum_M_over_D = X_mean * sum_m / denominator,
    n_sum_M = n * sum_m,
    V_sum_M = V * sum_m,
    X_sum_M = X_mean * sum_m
  )
  scaled <- ap_baseline_polynomial(n, V, X_mean, N / denominator)
  raw <- ap_baseline_polynomial(n, V, X_mean, N)
  names(raw) <- paste0("raw_", names(raw))
  cbind(base, scaled, raw)
}

ap_reduced_d_features <- function(X, component = c("full", "monly")) {
  component <- match.arg(component)
  shared <- c(
    "N", "n", "log_D", "sum_M", "sum_M2", "sum_M_sq_over_D",
    "sum_M_over_D", "n_sum_M_over_D", "n_sum_M",
    "base_0", "base_n", "base_n2",
    "raw_base_0", "raw_base_n", "raw_base_n2"
  )
  full_only <- c(
    "sum_Z", "sum_Z2", "sum_ZM", "sum_Z_sq_over_D",
    "sum_ZM_over_D", "sum_Z_over_D", "n_sum_Z_over_D"
  )
  keep <- if (component == "full") c(shared, full_only) else shared
  ap_assert(all(keep %in% names(X)),
            paste("Reduced", component, "density-ratio features are missing."))
  X[, keep, drop = FALSE]
}

ap_prepare_augmented_data <- function(train, arm, law, rho_full, rho_m,
                                      config, K, rep_id, fold_id) {
  arm_data <- train[train$A == arm, , drop = FALSE]
  ids <- sort(unique(arm_data$id))
  ap_assert(length(ids) >= 6L,
            "At least six training clusters per arm are required.")
  reference <- ap_fit_lambda_reference(
    arm_data,
    ap_hash_seed(config$base_seed, "lambda-fit", K, rep_id, fold_id, arm)
  )
  R <- config$reference_draws
  observed_full <- observed_m <- reference_full <- reference_m <- list()
  label <- numeric(0)
  classifier_id <- numeric(0)
  row_index <- 1L
  for (id in ids) {
    cluster <- arm_data[arm_data$id == id, , drop = FALSE]
    observed_full_row <- ap_d_full_features(cluster, cluster$M, rho_full)
    observed_m_row <- ap_d_m_features(cluster, cluster$M, rho_m)
    for (r in seq_len(R)) {
      reference_mediator <- ap_draw_lambda(
        reference, cluster, law,
        ap_hash_seed(config$base_seed, "lambda-draw", K, rep_id, fold_id,
                     arm, law, id, r)
      )
      observed_full[[row_index]] <- observed_full_row
      observed_m[[row_index]] <- observed_m_row
      reference_full[[row_index]] <- ap_d_full_features(
        cluster, reference_mediator, rho_full
      )
      reference_m[[row_index]] <- ap_d_m_features(
        cluster, reference_mediator, rho_m
      )
      label <- c(label, 0, 1)
      classifier_id <- c(classifier_id, id, id)
      row_index <- row_index + 1L
    }
  }
  list(
    label = label,
    id = classifier_id,
    full = rbind(do.call(rbind, observed_full),
                 do.call(rbind, reference_full)),
    monly = rbind(do.call(rbind, observed_m),
                  do.call(rbind, reference_m)),
    ordering = "blocked"
  )
}

ap_reorder_augmented <- function(prepared) {
  Rn <- nrow(prepared$full) / 2L
  index <- as.vector(rbind(seq_len(Rn), Rn + seq_len(Rn)))
  list(
    label = rep(c(0, 1), Rn),
    id = prepared$id,
    full = prepared$full[index, , drop = FALSE],
    monly = prepared$monly[index, , drop = FALSE]
  )
}

ap_grouped_foldid <- function(y, id, V = 5L, seed = 1L) {
  units <- unique(id)
  unit_y <- vapply(units, function(unit) mean(y[id == unit]), 0.0)
  strata <- ifelse(unit_y <= 0.25, "low",
                   ifelse(unit_y >= 0.75, "high", "mixed"))
  unit_fold <- integer(length(units))
  set.seed(seed)
  for (stratum in unique(strata)) {
    index <- which(strata == stratum)
    index <- sample(index, length(index))
    unit_fold[index] <- rep(seq_len(V), length.out = length(index))
  }
  setNames(unit_fold, units)[as.character(id)]
}

ap_classifier_matrix <- function(X) {
  X <- ap_model_frame(X)
  matrix <- as.matrix(X)
  storage.mode(matrix) <- "double"
  matrix
}

ap_fit_glmnet_classifier <- function(y, X, foldid, seed, alpha = 1) {
  set.seed(seed)
  fit <- glmnet::cv.glmnet(
    x = ap_classifier_matrix(X), y = as.numeric(y), family = "binomial",
    alpha = alpha, foldid = as.integer(foldid), type.measure = "deviance",
    nlambda = 40L, lambda.min.ratio = 0.02, standardize = TRUE,
    keep = TRUE, parallel = FALSE, control = list(maxit = 1000000L)
  )
  lambda_index <- which.min(abs(fit$lambda - fit$lambda.min))
  oof <- plogis(fit$fit.preval[, lambda_index])
  list(
    model = list(type = "ap_glmnet", fit = fit, names = names(X),
                 lambda = "lambda.min"),
    oof = as.numeric(oof)
  )
}

ap_fit_glm_classifier <- function(y, X, foldid) {
  folds <- sort(unique(foldid))
  oof <- rep(NA_real_, length(y))
  for (fold in folds) {
    validation <- foldid == fold
    fit <- ap_fit_base_model(
      y[!validation], X[!validation, , drop = FALSE], "binomial"
    )
    oof[validation] <- ap_predict_base_model(
      fit, X[validation, , drop = FALSE]
    )
  }
  list(model = ap_fit_base_model(y, X, "binomial"), oof = oof)
}

ap_fit_ranger_once <- function(y, X, seed) {
  data <- ap_model_frame(X)
  data$.outcome <- factor(as.integer(y), levels = c(0, 1))
  ranger::ranger(
    dependent.variable.name = ".outcome", data = data,
    probability = TRUE, num.trees = 300L,
    mtry = max(1L, floor(sqrt(ncol(X)))), min.node.size = 5L,
    sample.fraction = 0.8, replace = TRUE, importance = "none",
    num.threads = 1L, seed = seed, write.forest = TRUE
  )
}

ap_predict_ranger_once <- function(fit, X) {
  prediction <- predict(fit, data = ap_model_frame(X), num.threads = 1L)$predictions
  if (is.null(dim(prediction))) return(as.numeric(prediction))
  if ("1" %in% colnames(prediction)) return(as.numeric(prediction[, "1"]))
  as.numeric(prediction[, ncol(prediction)])
}

ap_fit_ranger_classifier <- function(y, X, foldid, seed) {
  folds <- sort(unique(foldid))
  oof <- rep(NA_real_, length(y))
  for (fold in folds) {
    validation <- foldid == fold
    fit_fold <- ap_fit_ranger_once(
      y[!validation], X[!validation, , drop = FALSE], seed + fold
    )
    oof[validation] <- ap_predict_ranger_once(
      fit_fold, X[validation, , drop = FALSE]
    )
  }
  final <- ap_fit_ranger_once(y, X, seed + 100L)
  list(
    model = list(type = "ap_ranger", fit = final, names = names(X)),
    oof = oof
  )
}

ap_fit_classifier_suite <- function(y, X, id, seed, linear_alpha = 1,
                                    linear_engine = "glmnet",
                                    sl_library = c("SL.glm", "SL.ranger"),
                                    requested = c(linear_engine, "ranger", "sl")) {
  X <- ap_model_frame(X)
  V <- max(2L, min(5L, length(unique(id))))
  foldid <- ap_grouped_foldid(y, id, V, seed)
  output <- list()
  if (linear_engine %in% requested) {
    linear <- if (linear_engine == "glm") {
      ap_fit_glm_classifier(y, X, foldid)
    } else {
      ap_fit_glmnet_classifier(
        y, X, foldid, seed + 10L, linear_alpha
      )
    }
    output[[linear_engine]] <- linear$model
  }
  if ("ranger" %in% requested) {
    output$ranger <- ap_fit_ranger_classifier(
      y, X, foldid, seed + 20L
    )$model
  }
  if ("sl" %in% requested) {
    ensemble <- ap_fit_model(
      y, X, "binomial", "ml", id, seed + 30L, sl_library
    )
    output$sl <- ensemble
    output$oof_weights <- if (ensemble$type == "SuperLearner") {
      stats::setNames(as.numeric(ensemble$fit$coef), names(ensemble$fit$coef))
    } else {
      stats::setNames(c(1, 0), sl_library)
    }
  }
  output
}

ap_predict_classifier <- function(model, X) {
  X <- ap_model_frame(X)
  missing <- setdiff(model$names, names(X))
  for (name in missing) X[[name]] <- 0
  X <- X[, model$names, drop = FALSE]
  if (model$type %in% c("base_lm", "base_glm")) {
    return(ap_predict_base_model(model, X))
  }
  if (model$type == "ap_glmnet") {
    return(as.numeric(predict(
      model$fit, newx = ap_classifier_matrix(X),
      s = model$lambda, type = "response"
    )))
  }
  if (model$type == "ap_ranger") {
    return(ap_predict_ranger_once(model$fit, X))
  }
  if (model$type == "SuperLearner") {
    return(ap_predict_model(model, X))
  }
  if (model$type == "ap_ensemble") {
    prediction <- Map(ap_predict_classifier, model$models, list(X))
    return(drop(cbind(prediction[[1L]], prediction[[2L]]) %*% model$weights))
  }
  ap_stop("Unknown classifier model type: ", model$type)
}

ap_fit_augmented_classifier <- function(prepared, classifier, seed,
                                        reduced = FALSE) {
  prepared <- ap_reorder_augmented(prepared)
  full <- if (reduced) ap_reduced_d_features(prepared$full, "full") else
    prepared$full
  monly <- if (reduced) ap_reduced_d_features(prepared$monly, "monly") else
    prepared$monly
  list(
    full = ap_fit_model(
      prepared$label, full, "binomial", classifier,
      prepared$id, seed
    ),
    monly = ap_fit_model(
      prepared$label, monly, "binomial", classifier,
      prepared$id, seed + 1L
    )
  )
}

ap_fit_augmented_suite <- function(prepared, seed, reduced = FALSE,
                                   sl_library = c("SL.glm", "SL.ranger"),
                                   classifier_names = c("glmnet", "ranger", "sl")) {
  prepared <- ap_reorder_augmented(prepared)
  full_X <- if (reduced) ap_reduced_d_features(prepared$full, "full") else
    prepared$full
  monly_X <- if (reduced) ap_reduced_d_features(prepared$monly, "monly") else
    prepared$monly
  full <- ap_fit_classifier_suite(
    prepared$label, full_X, prepared$id, seed, linear_alpha = 1,
    sl_library = sl_library, requested = classifier_names
  )
  monly <- ap_fit_classifier_suite(
    prepared$label, monly_X, prepared$id, seed + 1000L,
    linear_alpha = 1, sl_library = sl_library,
    requested = classifier_names
  )
  lapply(classifier_names, function(name) {
    list(full = full[[name]], monly = monly[[name]],
         weights = list(full = full$oof_weights, monly = monly$oof_weights))
  }) |>
    setNames(classifier_names)
}

ap_fit_d_models <- function(train, process, config, K, rep_id, fold_id) {
  models <- list(correct = list(), wrong = list())
  for (law in config$lambda_laws) {
    for (arm in 0:1) {
      prepared <- ap_prepare_augmented_data(
        train, arm, law, process$xi$rho, process$kappa$rho,
        config, K, rep_id, fold_id
      )
      for (status in c("correct", "wrong")) {
        suite <- ap_fit_augmented_suite(
          prepared,
          ap_hash_seed(config$base_seed, "d-suite", K, rep_id, fold_id,
                       law, arm, status),
          reduced = identical(status, "wrong"),
          sl_library = config$lowdim_s_library,
          classifier_names = config$classifier_names
        )
        for (classifier in config$classifier_names) {
          method <- paste("ld", law, classifier, sep = "_")
          if (is.null(models[[status]][[method]])) {
            models[[status]][[method]] <- list()
          }
          models[[status]][[method]][[as.character(arm)]] <-
            suite[[classifier]]
        }
      }
    }
  }
  models
}

ap_predict_d <- function(model, validation, rho_full, rho_m, config) {
  ids <- sort(unique(validation$id))
  full <- do.call(rbind, lapply(ids, function(id) {
    cluster <- validation[validation$id == id, , drop = FALSE]
    ap_d_full_features(cluster, cluster$M, rho_full)
  }))
  monly <- do.call(rbind, lapply(ids, function(id) {
    cluster <- validation[validation$id == id, , drop = FALSE]
    ap_d_m_features(cluster, cluster$M, rho_m)
  }))
  p_full <- ap_clip_prob(ap_predict_classifier(model$full, full), config)
  p_m <- ap_clip_prob(ap_predict_classifier(model$monly, monly), config)
  raw <- (p_full / (1 - p_full)) * ((1 - p_m) / p_m)
  list(
    value = setNames(ap_clip_ratio(raw, config), ids),
    raw = setNames(raw, ids),
    clipped = setNames(raw < config$ratio_lower | raw > config$ratio_upper, ids)
  )
}

ap_s_baseline_features <- function(n, V, X_mean, size_over_D) {
  data.frame(
    base_0 = size_over_D,
    base_n = size_over_D * n,
    base_V = size_over_D * V,
    base_X = size_over_D * X_mean
  )
}

ap_s_features_cluster <- function(data, rho_basis, reduced = FALSE) {
  cf <- ap_cluster_frame(data)
  M_sum <- as.numeric(rowsum(data$M, data$id, reorder = TRUE))
  N <- cf$N
  n <- ap_n_center(N)
  denominator <- 1 + (N - 1) * rho_basis
  base <- data.frame(
    M_sum_over_D = M_sum / denominator
  )
  out <- cbind(base, ap_s_baseline_features(
    n, cf$V, cf$X_mean, N / denominator
  ))
  if (reduced) out <- out[, c("M_sum_over_D", "base_0", "base_n"), drop = FALSE]
  out
}

ap_s_features_minus <- function(data, rho_basis, reduced = FALSE) {
  q <- data$N - 1
  n <- ap_n_center(data$N)
  M_sum <- q * data$M_peer
  denominator <- 1 + (q - 1) * rho_basis
  base <- data.frame(
    M_sum_over_D = M_sum / denominator
  )
  out <- cbind(base, ap_s_baseline_features(
    n, data$V, data$X_mean, q / denominator
  ))
  if (reduced) out <- out[, c("M_sum_over_D", "base_0", "base_n"), drop = FALSE]
  out
}

ap_fit_treatment_classifiers <- function(train, rho_basis, config, K,
                                         rep_id, fold_id) {
  cf <- ap_cluster_frame(train)
  output <- list()
  for (status in c("correct", "wrong")) {
    reduced <- identical(status, "wrong")
    cluster_suite <- ap_fit_classifier_suite(
      cf$A, ap_s_features_cluster(train, rho_basis, reduced), cf$id,
      ap_hash_seed(config$base_seed, "s-cluster", K, rep_id, fold_id, status),
      linear_engine = "glm", sl_library = config$lowdim_s_library,
      requested = "sl"
    )
    minus_suite <- ap_fit_classifier_suite(
      train$A, ap_s_features_minus(train, rho_basis, reduced), train$id,
      ap_hash_seed(config$base_seed, "s-minus", K, rep_id, fold_id, status),
      linear_engine = "glm", sl_library = config$lowdim_s_library,
      requested = "sl"
    )
    output[[status]] <- list(
      cluster = cluster_suite$sl, minus = minus_suite$sl,
      weights = list(cluster = cluster_suite$oof_weights,
                     minus = minus_suite$oof_weights)
    )
  }
  output
}

ap_predict_treatment_classifiers <- function(fitted, validation, rho_basis,
                                             config, reduced = FALSE) {
  cf <- ap_cluster_frame(validation)
  cluster_raw <- ap_predict_classifier(
    fitted$cluster, ap_s_features_cluster(validation, rho_basis, reduced)
  )
  minus_raw <- ap_predict_classifier(
    fitted$minus, ap_s_features_minus(validation, rho_basis, reduced)
  )
  list(
    cluster = setNames(ap_clip_prob(cluster_raw, config), cf$id),
    minus = ap_clip_prob(minus_raw, config),
    clipping = list(
      cluster = cluster_raw < config$prob_lower |
        cluster_raw > config$prob_upper,
      minus = minus_raw < config$prob_lower | minus_raw > config$prob_upper
    )
  )
}

ap_true_s_cluster <- function(data, kappa, config) {
  ids <- sort(unique(data$id))
  answer <- numeric(length(ids))
  pi <- config$dgp$treatment_prob
  for (i in seq_along(ids)) {
    cluster <- data[data$id == ids[[i]], , drop = FALSE]
    mean_1 <- ap_process_prediction(kappa, cluster, 1)
    mean_0 <- ap_process_prediction(kappa, cluster, 0)
    logit <- qlogis(pi) +
      ap_eq_log_density(cluster$M, mean_1, kappa$sigma, kappa$rho) -
      ap_eq_log_density(cluster$M, mean_0, kappa$sigma, kappa$rho)
    answer[[i]] <- plogis(logit)
  }
  setNames(ap_clip_prob(answer, config), ids)
}

ap_true_s_minus <- function(data, kappa, config) {
  pi <- config$dgp$treatment_prob
  answer <- numeric(nrow(data))
  for (id in sort(unique(data$id))) {
    index <- which(data$id == id)
    cluster <- data[index, , drop = FALSE]
    mean_1 <- ap_process_prediction(kappa, cluster, 1)
    mean_0 <- ap_process_prediction(kappa, cluster, 0)
    logit <- qlogis(pi) +
      ap_eq_leave_log_density(cluster$M, mean_1, kappa$sigma, kappa$rho) -
      ap_eq_leave_log_density(cluster$M, mean_0, kappa$sigma, kappa$rho)
    answer[index] <- plogis(logit)
  }
  ap_clip_prob(answer, config)
}

ap_true_d_cluster <- function(data, arm, kappa, xi, config) {
  ids <- sort(unique(data$id))
  value <- raw <- numeric(length(ids))
  for (i in seq_along(ids)) {
    cluster <- data[data$id == ids[[i]], , drop = FALSE]
    ratio <- ap_canonical_theta_ratio_cluster(
      cluster, kappa, xi, arm, arm, config
    )
    value[[i]] <- ratio$value
    raw[[i]] <- ratio$raw
  }
  list(
    value = setNames(value, ids), raw = setNames(raw, ids),
    clipped = setNames(raw < config$ratio_lower | raw > config$ratio_upper, ids)
  )
}

ap_kappa_scalar_ratio <- function(kappa, validation, config) {
  raw <- numeric(nrow(validation))
  for (id in sort(unique(validation$id))) {
    index <- which(validation$id == id)
    cluster <- validation[index, , drop = FALSE]
    mean_1 <- ap_process_prediction(kappa, cluster, 1)
    pieces <- ap_eq_marginal_conditional(
      cluster$M, mean_1, kappa$sigma, kappa$rho
    )
    raw[index] <- exp(ap_clip(
      pieces$marginal_log - pieces$conditional_log, -700, 700
    ))
  }
  list(
    value = ap_clip_ratio(raw, config), raw = raw,
    clipped = raw < config$ratio_lower | raw > config$ratio_upper
  )
}

ap_fit_lowdim_common <- function(train, config, K, rep_id, fold_id,
                                 fit_classifiers = TRUE) {
  eta <- list(correct = list(), wrong = list())
  for (correct in c(TRUE, FALSE)) {
    status <- ap_bool_name(correct)
    for (arm in 0:1) {
      eta[[status]][[as.character(arm)]] <- ap_fit_eta(
        train, arm, correct, "par",
        ap_hash_seed(config$base_seed, "lowdim", K, rep_id, fold_id,
                     "eta", status, arm)
      )
    }
  }
  phi <- ap_fit_process(
    train, "Z", TRUE, FALSE, "par",
    ap_hash_seed(config$base_seed, "lowdim", K, rep_id, fold_id, "phi")
  )
  xi <- ap_fit_process(
    train, "M", TRUE, TRUE, "par",
    ap_hash_seed(config$base_seed, "lowdim", K, rep_id, fold_id, "xi")
  )
  process <- list(
    phi = phi,
    xi = xi,
    kappa = ap_derive_compatible_kappa(phi, xi, train)
  )
  output <- list(eta = eta, process = process)
  if (fit_classifiers) {
    output$s <- ap_fit_treatment_classifiers(
      train, process$kappa$rho, config, K, rep_id, fold_id
    )
    output$d <- ap_fit_d_models(train, process, config, K, rep_id, fold_id)
  }
  output
}

ap_lowdim_validation_cache <- function(validation, fitted, config, K, rep_id,
                                       fold_id, oracle_ratios = FALSE) {
  ids <- sort(unique(validation$id))
  n <- nrow(validation)
  eta <- list(correct = list(), wrong = list())
  theta_parts <- tau_parts <- list()
  for (status in c("correct", "wrong")) {
    eta[[status]] <- list()
    for (arm in 0:1) {
      eta[[status]][[as.character(arm)]] <- ap_predict_eta(
        fitted$eta[[status]][[as.character(arm)]], validation
      )
    }
    theta_parts[[status]] <- list(
      theta11 = ap_partial_theta_par(
        fitted$eta[[status]][["1"]], validation,
        fitted$process$phi, fitted$process$kappa, 1, 1
      ),
      theta10 = ap_partial_theta_par(
        fitted$eta[[status]][["1"]], validation,
        fitted$process$phi, fitted$process$kappa, 1, 0
      ),
      theta00 = ap_partial_theta_par(
        fitted$eta[[status]][["0"]], validation,
        fitted$process$phi, fitted$process$kappa, 0, 0
      )
    )
    tau_parts[[status]] <- ap_partial_tau_par(
      fitted$eta[[status]][["1"]], validation,
      fitted$process$phi, fitted$process$kappa
    )
  }

  if (oracle_ratios) {
    xi <- ap_true_process_object("xi", config)
    d_true <- list(
      `0` = ap_true_d_cluster(
        validation, 0, fitted$process$kappa, xi, config
      ),
      `1` = ap_true_d_cluster(
        validation, 1, fitted$process$kappa, xi, config
      )
    )
    s_true <- list(
      cluster = ap_true_s_cluster(validation, fitted$process$kappa, config),
      minus = ap_true_s_minus(validation, fitted$process$kappa, config)
    )
    s_true$clipping <- list(
      cluster = s_true$cluster <= config$prob_lower |
        s_true$cluster >= config$prob_upper,
      minus = s_true$minus <= config$prob_lower |
        s_true$minus >= config$prob_upper
    )
    # Oracle mode is used only for exact score-identity checks.  Its wrong
    # branch is deliberately identical because misspecification is irrelevant
    # to the all-correct identity being tested.
    d_methods <- list(
      correct = list(oracle = d_true),
      wrong = list(oracle = d_true)
    )
    s <- list(correct = s_true, wrong = s_true)
  } else {
    s <- lapply(c("correct", "wrong"), function(status) {
      ap_predict_treatment_classifiers(
        fitted$s[[status]], validation, fitted$process$kappa$rho, config,
        reduced = identical(status, "wrong")
      )
    })
    names(s) <- c("correct", "wrong")
    d_methods <- lapply(c("correct", "wrong"), function(status) {
      lapply(fitted$d[[status]], function(method) {
        list(
          `0` = ap_predict_d(
            method[["0"]], validation, fitted$process$xi$rho,
            fitted$process$kappa$rho, config
          ),
          `1` = ap_predict_d(
            method[["1"]], validation, fitted$process$xi$rho,
            fitted$process$kappa$rho, config
          )
        )
      })
    })
    names(d_methods) <- c("correct", "wrong")
  }
  list(
    eta = eta, theta_parts = theta_parts, tau_parts = tau_parts,
    d = d_methods, s = s,
    kappa_scalar = ap_kappa_scalar_ratio(fitted$process$kappa, validation, config)
  )
}

ap_lowdim_theta_score <- function(validation, a, a_prime, eta, d_cluster,
                                  s_cluster, parts, config,
                                  stabilize = TRUE) {
  ids <- as.character(validation$id)
  d_row <- as.numeric(d_cluster[ids])
  s <- as.numeric(s_cluster[ids])
  s_num <- if (a_prime == 1) s else 1 - s
  s_den <- if (a == 1) s else 1 - s
  pi_a <- if (a == 1) config$dgp$treatment_prob else
    1 - config$dgp$treatment_prob
  pi_ap <- if (a_prime == 1) config$dgp$treatment_prob else
    1 - config$dgp$treatment_prob
  odds <- ap_clip_ratio(s_num / s_den, config)
  ratio <- ap_clip_ratio(d_row * (pi_a / pi_ap) * odds, config)
  ap_score_theta_fold(
    validation, a, a_prime, eta, ratio, parts, config, stabilize
  )
}

ap_lowdim_tau_score <- function(validation, eta, e_row, s_minus, parts,
                                config, stabilize = TRUE) {
  pi <- config$dgp$treatment_prob
  odds <- ap_clip_ratio((1 - s_minus) / s_minus, config)
  ratio <- ap_clip_ratio(e_row * (pi / (1 - pi)) * odds, config)
  ap_score_tau_fold(validation, eta, ratio, parts, config, stabilize)
}

ap_score_lowdim_cache <- function(validation, cache, scenario, d_method,
                                  config, stabilize = TRUE) {
  theta_cfg <- scenario$theta
  tau_cfg <- scenario$tau
  theta_status <- ap_bool_name(theta_cfg$eta)
  tau_status <- ap_bool_name(tau_cfg$eta)
  d_theta <- cache$d[[ap_bool_name(theta_cfg$d)]][[d_method]]
  d_tau <- cache$d[[ap_bool_name(tau_cfg$e)]][[d_method]]
  s_cluster <- cache$s[[ap_bool_name(theta_cfg$s)]]$cluster
  s_minus <- cache$s[[ap_bool_name(tau_cfg$s)]]$minus
  definitions <- list(
    theta11 = c(a = 1, a_prime = 1),
    theta10 = c(a = 1, a_prime = 0),
    theta00 = c(a = 0, a_prime = 0)
  )
  components <- list()
  for (component in names(definitions)) {
    definition <- definitions[[component]]
    d_value <- d_theta[[as.character(definition[["a"]])]]$value
    components[[component]] <- ap_lowdim_theta_score(
      validation, definition[["a"]], definition[["a_prime"]],
      cache$eta[[theta_status]][[as.character(definition[["a"]])]],
      d_value, s_cluster, cache$theta_parts[[theta_status]][[component]],
      config, stabilize
    )
  }
  d1_row <- as.numeric(d_tau[["1"]]$value[as.character(validation$id)])
  e_row <- ap_clip_ratio(d1_row * cache$kappa_scalar$value, config)
  components$tau <- ap_lowdim_tau_score(
    validation, cache$eta[[tau_status]][["1"]], e_row, s_minus,
    cache$tau_parts[[tau_status]], config, stabilize
  )
  components
}

ap_lowdim_fold_diagnostics <- function(validation, cache, method, scenario,
                                       config, K, rep_id, fold_id) {
  rows <- list()
  index <- 1L
  add <- function(component, quantity, detail, indicator) {
    indicator <- as.logical(indicator)
    rows[[index]] <<- data.frame(
      study = "lowdim", K = K, rep = rep_id,
      method = paste0(method, "_stabilized"), fold = fold_id,
      scenario = scenario, component = component, quantity = quantity,
      detail = detail, n = length(indicator), truncated = sum(indicator),
      truncated_fraction = mean(indicator)
    )
    index <<- index + 1L
  }
  scenario_cfg <- ap_lowdim_scenarios()[[scenario]]
  ids <- as.character(validation$id)
  theta_d_status <- ap_bool_name(scenario_cfg$theta$d)
  tau_d_status <- ap_bool_name(scenario_cfg$tau$e)
  theta_s_status <- ap_bool_name(scenario_cfg$theta$s)
  tau_s_status <- ap_bool_name(scenario_cfg$tau$s)
  s_cluster <- cache$s[[theta_s_status]]$cluster
  s_minus <- cache$s[[tau_s_status]]$minus
  add("s", "classifier_probability", paste0("cluster_", theta_s_status),
      cache$s[[theta_s_status]]$clipping$cluster)
  add("s_minus", "classifier_probability", paste0("leave_one_", tau_s_status),
      cache$s[[tau_s_status]]$clipping$minus)
  add("kappa_j_over_kappa_star", "density_ratio", "shared",
      cache$kappa_scalar$clipped)
  definitions <- list(
    theta11 = c(a = 1, a_prime = 1),
    theta10 = c(a = 1, a_prime = 0),
    theta00 = c(a = 0, a_prime = 0)
  )
  d_theta <- cache$d[[theta_d_status]][[method]]
  for (component in names(definitions)) {
    definition <- definitions[[component]]
    arm <- as.character(definition[["a"]])
    d_value <- d_theta[[arm]]$value
    add(component, "density_ratio",
        paste0("d_", theta_d_status, "_arm", arm), d_theta[[arm]]$clipped)
    s_num <- if (definition[["a_prime"]] == 1) s_cluster else 1 - s_cluster
    s_den <- if (definition[["a"]] == 1) s_cluster else 1 - s_cluster
    odds <- ap_clip_ratio(s_num / s_den, config)
    raw_ratio <- as.numeric(d_value) * odds
    add(component, "score_ratio", "theta",
        raw_ratio < config$ratio_lower | raw_ratio > config$ratio_upper)
  }
  d_tau <- cache$d[[tau_d_status]][[method]]
  d1_row <- as.numeric(d_tau[["1"]]$value[ids])
  e_raw <- d1_row * cache$kappa_scalar$value
  e <- ap_clip_ratio(e_raw, config)
  add("tau", "density_ratio", paste0("e_", tau_d_status),
      e_raw < config$ratio_lower | e_raw > config$ratio_upper)
  odds_tau <- ap_clip_ratio((1 - s_minus) / s_minus, config)
  raw_tau <- e * odds_tau
  add("tau", "score_ratio", "tau",
      raw_tau < config$ratio_lower | raw_tau > config$ratio_upper)
  do.call(rbind, rows)
}

ap_run_lowdim_eif <- function(data, K, rep_id, config, stabilize = TRUE) {
  scenarios <- ap_lowdim_scenarios()
  methods <- as.vector(outer(
    paste0("ld_", config$lambda_laws), config$classifier_names,
    function(x, y) paste(x, y, sep = "_")
  ))
  ids <- sort(unique(data$id))
  components <- lapply(methods, function(method) {
    lapply(scenarios, function(x) list(
      theta11 = ap_empty_component(ids), theta10 = ap_empty_component(ids),
      theta00 = ap_empty_component(ids), tau = ap_empty_component(ids)
    ))
  })
  names(components) <- methods
  diagnostics <- list()
  diagnostic_index <- 1L
  fold <- ap_make_folds(
    data, config$folds,
    ap_hash_seed(config$base_seed, "fold", K, rep_id)
  )
  for (fold_id in seq_len(config$folds)) {
    validation_ids <- as.numeric(names(fold)[fold == fold_id])
    validation <- data[data$id %in% validation_ids, , drop = FALSE]
    train <- data[!data$id %in% validation_ids, , drop = FALSE]
    fitted <- ap_fit_lowdim_common(train, config, K, rep_id, fold_id, TRUE)
    cache <- ap_lowdim_validation_cache(
      validation, fitted, config, K, rep_id, fold_id, FALSE
    )
    for (method in methods) {
      for (scenario_id in names(scenarios)) {
        scored <- ap_score_lowdim_cache(
          validation, cache, scenarios[[scenario_id]], method,
          config, stabilize
        )
        components[[method]][[scenario_id]] <- ap_merge_fold_component(
          components[[method]][[scenario_id]], scored
        )
        diagnostics[[diagnostic_index]] <- ap_lowdim_fold_diagnostics(
          validation, cache, method, scenario_id, config, K, rep_id, fold_id
        )
        diagnostic_index <- diagnostic_index + 1L
      }
    }
  }
  cf <- ap_cluster_frame(data)
  N_by_id <- setNames(cf$N, cf$id)
  rows <- list()
  index <- 1L
  for (method in methods) {
    for (scenario_id in names(scenarios)) {
      effect_rows <- ap_effect_rows(
        components[[method]][[scenario_id]], N_by_id, config
      )
      rows[[index]] <- ap_finalize_rows(
        effect_rows, "lowdim", K, rep_id, scenario_id,
        paste0(method, if (stabilize) "_stabilized" else ""),
        scenarios[[scenario_id]]
      )
      index <- index + 1L
    }
  }
  list(rows = do.call(rbind, rows), components = components,
       diagnostics = do.call(rbind, diagnostics), fold = fold)
}

ap_score_lowdim_oracle_prefit <- function(validation, wrong_eta_fit, config,
                                          stabilize = TRUE,
                                          wrong_d = NULL, wrong_s = NULL) {
  fitted <- list(
    eta = list(
      correct = list(
        `0` = ap_true_eta_object(0, config),
        `1` = ap_true_eta_object(1, config)
      ),
      wrong = wrong_eta_fit
    ),
    process = list(
      phi = ap_true_process_object("phi", config),
      kappa = ap_true_process_object("kappa", config)
    )
  )
  cache <- ap_lowdim_validation_cache(
    validation, fitted, config, nrow(ap_cluster_frame(validation)),
    0L, 0L, TRUE
  )
  if (!is.null(wrong_d)) cache$d$wrong$oracle <- wrong_d
  if (!is.null(wrong_s)) cache$s$wrong <- wrong_s
  cf <- ap_cluster_frame(validation)
  N_by_id <- setNames(cf$N, cf$id)
  output <- list()
  for (scenario in ap_lowdim_scenarios()) {
    components <- ap_score_lowdim_cache(
      validation, cache, scenario, "oracle", config, stabilize
    )
    output[[scenario$id]] <- list(
      effects = ap_effect_rows(components, N_by_id, config),
      components = components
    )
  }
  output
}
