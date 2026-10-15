#!/usr/bin/env Rscript

source("annals_core.R")
source("annals_canonical.R")
source("annals_lowdim.R")

config <- ap_default_config()
ap_validate_config(config)

cat("[1] Fixed truncation constants\n")
ap_assert(config$prob_lower == 0.025 && config$prob_upper == 0.975,
          "Probability truncation changed.")
ap_assert(config$ratio_lower == 0.05 && config$ratio_upper == 20,
          "Ratio truncation changed.")
ap_assert(ap_ml_integration_draws(config, 50L) == config$mc_draws_ml,
          "The ML draw lower bound is not respected.")
ap_assert(ap_ml_integration_draws(config, 150L) == 150L,
          "ML integration draws must increase with K.")

cat("[1b] Fixed SL.glm + SL.ranger libraries\n")
ap_assert(identical(config$canonical_ml_library, c("SL.glm", "SL.ranger")),
          "Canonical ML library is not SL.glm + SL.ranger.")
ap_assert(identical(config$lowdim_s_library, c("SL.glm", "SL.ranger")),
          "Low-dimensional SL library is not SL.glm + SL.ranger.")
set.seed(1000)
smoke_id <- rep(seq_len(30L), each = 2L)
smoke_y <- rep(rep(c(0, 1), each = 15L), each = 2L)
smoke_X <- data.frame(x1 = rnorm(60L), x2 = rnorm(60L))
smoke_suite <- ap_fit_classifier_suite(
  smoke_y, smoke_X, smoke_id, 1000L, linear_engine = "glm",
  sl_library = config$lowdim_s_library
)
ap_assert(smoke_suite$sl$type == "SuperLearner" &&
            identical(smoke_suite$sl$library, c("SL.glm", "SL.ranger")),
          "Low-dimensional SL fit did not use the required two-wrapper library.")

cat("[1c] Prespecified vanishing-ridge separation fallback\n")
separation_id <- rep(seq_len(40L), each = 2L)
separation_x <- rep(rep(c(-1, 1), each = 20L), each = 2L)
separation_y <- as.numeric(separation_x > 0)
separation_X <- data.frame(
  x = separation_x,
  z = seq_along(separation_x) / length(separation_x)
)
separation_fit <- ap_fit_model(
  separation_y, separation_X, "binomial", "ml", separation_id, 1001L,
  c("SL.glm", "SL.ranger")
)
separation_prediction <- ap_predict_model(separation_fit, separation_X)
ap_assert(separation_fit$type == "SuperLearner" &&
            separation_fit$ridge_fallback_count > 0L,
          "Complete separation did not trigger the ridge fallback.")
ap_assert(identical(separation_fit$ridge_lambda, 1 / 40),
          "The fallback penalty is not one over the training-cluster count.")
ap_assert(all(is.finite(separation_prediction)),
          "The ridge fallback produced non-finite predictions.")
regular_id <- rep(seq_len(50L), each = 2L)
regular_y <- rep(c(0, 1), 50L)
regular_index <- seq_along(regular_y)
regular_X <- data.frame(
  x = sin(regular_index / 7),
  z = cos(regular_index / 11)
)
regular_fit <- ap_fit_model(
  regular_y, regular_X, "binomial", "ml", regular_id, 1002L,
  c("SL.glm", "SL.ranger")
)
ap_assert(identical(regular_fit$ridge_fallback_count, 0L),
          "A regular converged logistic GLM unnecessarily used ridge.")

ml_smoke_data <- ap_generate_data(30L, config, 1000L)
ml_smoke <- ap_fit_canonical_nuisances(
  ml_smoke_data, "ml", config, 30L, 1000L, 0L
)
ml_models <- c(
  unlist(lapply(ml_smoke$eta, unname), recursive = FALSE),
  unname(ml_smoke$phi), unname(ml_smoke$kappa), unname(ml_smoke$xi)
)
ap_assert(all(vapply(ml_models, function(fit) {
  fit$model$type == "SuperLearner" &&
    identical(fit$model$library, c("SL.glm", "SL.ranger"))
}, logical(1))), "A canonical ML nuisance omitted a required learner.")
ap_assert(!identical(ml_smoke$kappa$correct$type, "derived_kappa"),
          "Canonical ML kappa must be fitted directly with Super Learner.")

cat("[2] Deterministic truth quadrature\n")
truth_15 <- ap_truth(config, 15L)
truth_25 <- ap_truth(config, 25L)
ap_assert(max(abs(truth_15 - truth_25)) < 1e-10,
          "Truth quadrature is not numerically exact.")

cat("[3] Exchangeable Gaussian density formulas\n")
set.seed(1001)
for (n in c(3L, 10L, 20L)) {
  mean <- rnorm(n)
  y <- rnorm(n)
  sigma <- 1.3
  rho <- 0.08
  covariance <- sigma^2 * ((1 - rho) * diag(n) + rho * matrix(1, n, n))
  direct <- -n * log(2 * pi) / 2 - determinant(covariance, logarithm = TRUE)$modulus / 2 -
    drop(crossprod(y - mean, solve(covariance, y - mean))) / 2
  formula <- ap_eq_log_density(y, mean, sigma, rho)
  ap_assert(abs(direct - formula) < 1e-10,
            "Full exchangeable log density is incorrect.")
  leave <- ap_eq_leave_log_density(y, mean, sigma, rho)
  for (j in seq_len(n)) {
    keep <- setdiff(seq_len(n), j)
    direct_leave <- -length(keep) * log(2 * pi) / 2 -
      determinant(covariance[keep, keep], logarithm = TRUE)$modulus / 2 -
      drop(crossprod(y[keep] - mean[keep],
        solve(covariance[keep, keep], y[keep] - mean[keep]))) / 2
    ap_assert(abs(direct_leave - leave[[j]]) < 1e-10,
              "Leave-one exchangeable log density is incorrect.")
  }
}

cat("[4] Structural Gaussian compatibility\n")
d <- config$dgp
kappa <- ap_true_process_object("kappa", config)
expected_sigma <- sqrt(d$gamma_zm^2 * d$sigma_z^2 + d$sigma_m_error^2)
ap_assert(abs(kappa$sigma - expected_sigma) < 1e-14,
          "Marginal M variance is incompatible with M given Z.")
ap_assert(abs(kappa$rho - d$process_rho) < 1e-14,
          "Marginal M correlation is incompatible with M given Z.")

cat("[5] Fitted Gaussian-law compatibility\n")
compatibility_data <- ap_generate_data(500L, config, 1007L)
compatible <- ap_fit_canonical_nuisances(
  compatibility_data, "par", config, 500L, 1007L, 0L
)
for (arm in 0:1) {
  marginal_mean <- ap_process_prediction(
    compatible$kappa$correct, compatibility_data, arm
  )
  at_mean <- compatibility_data
  at_mean$Z <- ap_process_prediction(
    compatible$phi$correct, compatibility_data, arm
  )
  integrated_conditional_mean <- ap_process_prediction(
    compatible$xi$correct, at_mean, arm
  )
  ap_assert(max(abs(marginal_mean - integrated_conditional_mean)) < 1e-12,
            "Fitted marginal and conditional mediator means are incompatible.")
}
implied_variance <- compatible$kappa$correct$gamma^2 *
  compatible$phi$correct$sigma^2 + compatible$xi$correct$sigma^2
ap_assert(abs(compatible$kappa$correct$sigma^2 - implied_variance) < 1e-12,
          "Fitted marginal mediator variance is incompatible.")

cat("[6] Fitted augmented-classification feature spans\n")
for (arm in 0:1) {
  arm_data <- compatibility_data[compatibility_data$A == arm, , drop = FALSE]
  reference <- ap_fit_lambda_reference(arm_data, 1100L + arm)
  ids <- unique(arm_data$id)
  full_X <- monly_X <- vector("list", length(ids))
  full_log_odds <- monly_log_odds <- numeric(length(ids))
  for (i in seq_along(ids)) {
    cluster <- arm_data[arm_data$id == ids[[i]], , drop = FALSE]
    lambda_mean <- ap_predict_model(reference$model, ap_lambda_features(cluster))
    log_lambda <- sum(dnorm(
      cluster$M, lambda_mean, reference$sigma, log = TRUE
    ))
    full_log_odds[[i]] <- log_lambda - ap_eq_log_density(
      cluster$M,
      ap_process_prediction(compatible$xi$correct, cluster, arm),
      compatible$xi$correct$sigma, compatible$xi$correct$rho
    )
    monly_log_odds[[i]] <- log_lambda - ap_eq_log_density(
      cluster$M,
      ap_process_prediction(compatible$kappa$correct, cluster, arm),
      compatible$kappa$correct$sigma, compatible$kappa$correct$rho
    )
    full_X[[i]] <- ap_d_full_features(
      cluster, cluster$M, compatible$xi$correct$rho
    )
    monly_X[[i]] <- ap_d_m_features(
      cluster, cluster$M, compatible$kappa$correct$rho
    )
  }
  full_X <- do.call(rbind, full_X)
  monly_X <- do.call(rbind, monly_X)
  full_fit <- ap_fit_base_model(full_log_odds, full_X)
  monly_fit <- ap_fit_base_model(monly_log_odds, monly_X)
  full_rmse <- sqrt(mean((
    ap_predict_base_model(full_fit, full_X) - full_log_odds
  )^2))
  monly_rmse <- sqrt(mean((
    ap_predict_base_model(monly_fit, monly_X) - monly_log_odds
  )^2))
  ap_assert(max(full_rmse, monly_rmse) < 1e-8,
            "Fitted augmented-classification feature span is incorrect.")
}

cat("[7] Noiseless outcome-model span\n")
data <- ap_generate_data(500L, config, 1002L, include_y = FALSE)
for (arm in 0:1) {
  fit <- ap_fit_eta(data, arm, TRUE, "par", 1L)
  index <- data$A == arm
  error <- ap_predict_eta(fit, data[index, , drop = FALSE]) - data$eta0[index]
  ap_assert(max(abs(error)) < 1e-9,
            "Correct outcome regression does not span eta0.")
}

cat("[8] Treatment-stratified cluster folds\n")
fold <- ap_make_folds(data, config$folds, 1003L)
ap_assert(length(fold) == length(unique(data$id)), "Fold map is not cluster-level.")
ap_assert(all(table(fold) > 0), "At least one validation fold is empty.")
for (a in 0:1) {
  arm_ids <- ap_cluster_frame(data)$id[ap_cluster_frame(data)$A == a]
  ap_assert(length(unique(fold[as.character(arm_ids)])) == config$folds,
            "Treatment-stratified folds omit an arm in a fold.")
}

cat("[9] Stabilization estimating equation\n")
set.seed(1004)
Y <- rnorm(200)
eta <- rnorm(200)
clever <- runif(200, 0.05, 2)
epsilon <- ap_stabilization_epsilon(Y, eta, clever)
ap_assert(abs(sum(clever * (Y - eta - epsilon))) < 1e-10,
          "Stabilization fluctuation does not solve its estimating equation.")

cat("[10] Component-specific low-dimensional theory flags\n")
scenarios <- ap_lowdim_scenarios()
ap_assert(!scenarios[["LD-T4"]]$point_consistency_expected[["IDE"]],
          "LD-T4 IDE should not be protected.")
ap_assert(!scenarios[["LD-T4"]]$point_consistency_expected[["IIME"]],
          "LD-T4 IIME should not be protected.")
ap_assert(scenarios[["LD-T4"]]$point_consistency_expected[["ISME"]],
          "LD-T4 ISME is protected because s(a)/s(a)=1.")

cat("[11] True-nuisance score identity\n")
identity_data <- ap_generate_data(150L, config, 1005L, include_y = FALSE)
canonical_fit <- ap_fit_canonical_nuisances(
  identity_data, "par", config, 150L, 1005L, 0L
)
canonical_fit$eta$correct <- list(
  `0` = ap_true_eta_object(0, config),
  `1` = ap_true_eta_object(1, config)
)
canonical_fit$phi$correct <- ap_true_process_object("phi", config)
canonical_fit$kappa$correct <- ap_true_process_object("kappa", config)
canonical_fit$xi$correct <- ap_true_process_object("xi", config)
canonical_cache <- ap_canonical_validation_cache(
  identity_data, canonical_fit, "par", config, 150L, 1005L, 0L
)
canonical_score <- ap_score_canonical_cache(
  identity_data, canonical_cache, ap_canonical_scenarios()[["C-1111"]],
  config, stabilize = FALSE
)
wrong_eta <- list(
  `0` = ap_fit_eta(identity_data, 0, FALSE),
  `1` = ap_fit_eta(identity_data, 1, FALSE)
)
lowdim_fit <- list(
  eta = list(
    correct = list(
      `0` = ap_true_eta_object(0, config),
      `1` = ap_true_eta_object(1, config)
    ),
    wrong = wrong_eta
  ),
  process = list(
    phi = ap_true_process_object("phi", config),
    kappa = ap_true_process_object("kappa", config)
  )
)
lowdim_cache <- ap_lowdim_validation_cache(
  identity_data, lowdim_fit, config, 150L, 1005L, 0L,
  oracle_ratios = TRUE
)
lowdim_score <- ap_score_lowdim_cache(
  identity_data, lowdim_cache, scenarios[["LD-AC"]], "oracle",
  config, stabilize = FALSE
)
for (component in names(canonical_score)) {
  ap_assert(max(abs(canonical_score[[component]]$C -
    lowdim_score[[component]]$C)) < 1e-10,
  paste("Cluster score identity failed for", component))
  ap_assert(max(abs(canonical_score[[component]]$I -
    lowdim_score[[component]]$I)) < 1e-10,
  paste("Individual score identity failed for", component))
}

cat("[12] Reproducible data generation\n")
first <- ap_generate_data(20L, config, 1006L)
second <- ap_generate_data(20L, config, 1006L)
ap_assert(identical(first, second), "DGP seeds are not reproducible.")

cat("[13] Reduced outcome model remains treatment-contrast nondegenerate\n")
reduced_eta_0 <- ap_fit_eta(data, 0, FALSE, "par", 1201L)
reduced_eta_1 <- ap_fit_eta(data, 1, FALSE, "par", 1202L)
ap_assert(
  max(abs(ap_predict_eta(reduced_eta_1, data) -
          ap_predict_eta(reduced_eta_0, data))) > 1e-3,
  "Reduced arm-specific outcome regressions collapse to one common fit."
)
gcomp_matrix <- ap_gcomp_point_matrix(data, config, 1203L)
ap_assert(any(abs(gcomp_matrix[substr(rownames(gcomp_matrix), 1L, 1L) == "0", ]) >
                1e-3),
          "Reduced eta forces all g-computation contrasts to zero.")

cat("[14] Reduced process laws retain defining conditioning arguments\n")
xi_reduced_names <- names(ap_process_features(
  data, correct = FALSE, conditional = TRUE
))
ap_assert(all(c("Z", "A", "n") %in% xi_reduced_names),
          "Reduced xi must retain Z, treatment, and cluster size.")
ap_assert(!any(c("V", "X_mean", "A_n") %in% xi_reduced_names),
          "Reduced xi unexpectedly retains omitted baseline terms.")

cat("[15] Reduced classifiers retain core variables but omit baselines\n")
cluster <- data[data$id == data$id[[1L]], , drop = FALSE]
d_full_reduced <- ap_reduced_d_features(
  ap_d_full_features(cluster, cluster$M, kappa$rho), "full"
)
d_m_reduced <- ap_reduced_d_features(
  ap_d_m_features(cluster, cluster$M, kappa$rho), "monly"
)
ap_assert(all(c("N", "sum_Z", "sum_M", "sum_ZM") %in%
                names(d_full_reduced)),
          "Reduced full classifier omits a defining process summary.")
ap_assert(all(c("N", "sum_M", "sum_M2") %in% names(d_m_reduced)),
          "Reduced mediator classifier omits a defining mediator summary.")
ap_assert(!any(grepl("(^|_)V|(^|_)X", names(d_full_reduced))),
          "Reduced full classifier retains an omitted baseline feature.")
s_reduced_names <- names(ap_s_features_cluster(data, kappa$rho, TRUE))
ap_assert(identical(s_reduced_names,
                    c("M_sum_over_D", "base_0", "base_n")),
          "Reduced treatment classifier has the wrong feature set.")

cat("All implementation checks passed.\n")
