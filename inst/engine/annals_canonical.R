# Canonical EIF estimators and parametric g-computation.

ap_empty_component <- function(ids) {
  list(
    C = setNames(rep(NA_real_, length(ids)), as.character(ids)),
    I = setNames(rep(NA_real_, length(ids)), as.character(ids))
  )
}

ap_bool_name <- function(x) if (isTRUE(x)) "correct" else "wrong"

ap_fit_canonical_nuisances <- function(train, learner, config, K, rep_id,
                                       fold_id = 0L) {
  eta <- list(correct = list(), wrong = list())
  for (correct in c(TRUE, FALSE)) {
    status <- ap_bool_name(correct)
    for (arm in 0:1) {
      eta[[status]][[as.character(arm)]] <- ap_fit_eta(
        train, arm, correct, learner,
        ap_hash_seed(config$base_seed, "canonical", learner, K, rep_id,
                     fold_id, "eta", status, arm),
        config$canonical_ml_library
      )
    }
  }
  fit_process_block <- function(outcome, conditional, correct) {
    status <- ap_bool_name(correct)
    ap_fit_process(
      train, outcome, correct, conditional, learner,
      ap_hash_seed(config$base_seed, "canonical", learner, K, rep_id,
                   fold_id, outcome, conditional, status),
      config$canonical_ml_library
    )
  }
  phi_correct <- fit_process_block("Z", FALSE, TRUE)
  phi_wrong <- fit_process_block("Z", FALSE, FALSE)
  xi_correct <- fit_process_block("M", TRUE, TRUE)
  xi_wrong <- fit_process_block("M", TRUE, FALSE)
  # The parametric fit retains the exact compatible Gaussian construction.
  # With SL.ranger in the ML library, xi need not be linear in Z, so its
  # marginal mediator law cannot be obtained by the linear-Gaussian identity.
  # Canonical theory permits kappa to be fitted directly; use the same ML
  # library on its correctly specified marginal feature set.
  kappa_correct <- if (learner == "par") {
    ap_derive_compatible_kappa(phi_correct, xi_correct, train)
  } else {
    fit_process_block("M", FALSE, TRUE)
  }
  kappa_wrong <- fit_process_block("M", FALSE, FALSE)
  list(
    eta = eta,
    phi = list(
      correct = phi_correct,
      wrong = phi_wrong
    ),
    kappa = list(
      correct = kappa_correct,
      wrong = kappa_wrong
    ),
    xi = list(
      correct = xi_correct,
      wrong = xi_wrong
    )
  )
}

ap_canonical_theta_ratio_cluster <- function(cluster, kappa, xi, a, a_prime,
                                              config) {
  mean_num <- ap_process_prediction(kappa, cluster, a_prime)
  mean_den <- ap_process_prediction(xi, cluster, a)
  log_ratio <- ap_eq_log_density(cluster$M, mean_num, kappa$sigma, kappa$rho) -
    ap_eq_log_density(cluster$M, mean_den, xi$sigma, xi$rho)
  raw <- exp(ap_clip(log_ratio, -700, 700))
  list(value = ap_clip_ratio(raw, config), raw = raw,
       clipped = raw < config$ratio_lower || raw > config$ratio_upper)
}

ap_canonical_tau_ratio_cluster <- function(cluster, kappa, xi, config) {
  mean_1 <- ap_process_prediction(kappa, cluster, 1)
  mean_0 <- ap_process_prediction(kappa, cluster, 0)
  mean_xi <- ap_process_prediction(xi, cluster, 1)
  log_marginal_1 <- dnorm(cluster$M, mean_1, kappa$sigma, log = TRUE)
  log_leave_0 <- ap_eq_leave_log_density(
    cluster$M, mean_0, kappa$sigma, kappa$rho
  )
  log_denominator <- ap_eq_log_density(
    cluster$M, mean_xi, xi$sigma, xi$rho
  )
  raw <- exp(ap_clip(log_marginal_1 + log_leave_0 - log_denominator, -700, 700))
  list(value = ap_clip_ratio(raw, config), raw = raw,
       clipped = raw < config$ratio_lower | raw > config$ratio_upper)
}

ap_draw_process <- function(process, cluster, arm, B, seed) {
  N <- nrow(cluster)
  mean <- ap_process_prediction(process, cluster, arm)
  covariance <- process$sigma^2 *
    ((1 - process$rho) * diag(N) + process$rho * matrix(1, N, N))
  eigenvalue <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalue) <= 1e-10) {
    covariance <- covariance + diag(abs(min(eigenvalue)) + 1e-8, N)
  }
  set.seed(seed)
  error <- matrix(rnorm(B * N), B, N) %*% chol(covariance)
  sweep(error, 2L, mean, `+`)
}

ap_stack_cluster <- function(cluster, B) {
  cluster[rep(seq_len(nrow(cluster)), times = B), , drop = FALSE]
}

ap_flatten_draws <- function(draws) as.vector(t(draws))

ap_peer_draws <- function(draws) {
  (rowSums(draws) - draws) / (ncol(draws) - 1)
}

ap_average_ml_prediction <- function(eta_fit, newdata, B, N) {
  prediction <- ap_predict_eta(eta_fit, newdata)
  colMeans(matrix(prediction, nrow = B, ncol = N, byrow = TRUE))
}

ap_partial_theta_ml <- function(eta_fit, cluster, M_draw, Z_draw) {
  B <- nrow(M_draw)
  N <- nrow(cluster)
  stacked <- ap_stack_cluster(cluster, B)
  observed_M <- rep(cluster$M, times = B)
  observed_M_peer <- rep(cluster$M_peer, times = B)
  observed_Z <- rep(cluster$Z, times = B)
  observed_Z_peer <- rep(cluster$Z_peer, times = B)

  u1_data <- stacked
  u1_data$M <- ap_flatten_draws(M_draw)
  u1_data$M_peer <- ap_flatten_draws(ap_peer_draws(M_draw))
  u1_data$Z <- observed_Z
  u1_data$Z_peer <- observed_Z_peer

  u2_data <- stacked
  u2_data$M <- observed_M
  u2_data$M_peer <- observed_M_peer
  u2_data$Z <- ap_flatten_draws(Z_draw)
  u2_data$Z_peer <- ap_flatten_draws(ap_peer_draws(Z_draw))

  u3_data <- u1_data
  u3_data$Z <- ap_flatten_draws(Z_draw)
  u3_data$Z_peer <- ap_flatten_draws(ap_peer_draws(Z_draw))

  list(
    u1 = ap_average_ml_prediction(eta_fit, u1_data, B, N),
    u2 = ap_average_ml_prediction(eta_fit, u2_data, B, N),
    u3 = ap_average_ml_prediction(eta_fit, u3_data, B, N)
  )
}

ap_partial_tau_ml <- function(eta_fit, cluster, M0, M1, Z1) {
  B <- nrow(M0)
  N <- nrow(cluster)
  stacked <- ap_stack_cluster(cluster, B)
  observed_M <- rep(cluster$M, times = B)
  observed_M_peer <- rep(cluster$M_peer, times = B)
  observed_Z <- rep(cluster$Z, times = B)
  observed_Z_peer <- rep(cluster$Z_peer, times = B)
  own_M1 <- ap_flatten_draws(M1)
  peer_M0 <- ap_flatten_draws(ap_peer_draws(M0))
  draw_Z1 <- ap_flatten_draws(Z1)
  peer_Z1 <- ap_flatten_draws(ap_peer_draws(Z1))

  r1_data <- stacked
  r1_data$M <- own_M1
  r1_data$M_peer <- peer_M0
  r1_data$Z <- observed_Z
  r1_data$Z_peer <- observed_Z_peer

  r2_data <- stacked
  r2_data$M <- observed_M
  r2_data$M_peer <- peer_M0
  r2_data$Z <- draw_Z1
  r2_data$Z_peer <- peer_Z1

  r3_data <- stacked
  r3_data$M <- own_M1
  r3_data$M_peer <- observed_M_peer
  r3_data$Z <- draw_Z1
  r3_data$Z_peer <- peer_Z1

  r4_data <- stacked
  r4_data$M <- own_M1
  r4_data$M_peer <- peer_M0
  r4_data$Z <- draw_Z1
  r4_data$Z_peer <- peer_Z1

  list(
    r1 = ap_average_ml_prediction(eta_fit, r1_data, B, N),
    r2 = ap_average_ml_prediction(eta_fit, r2_data, B, N),
    r3 = ap_average_ml_prediction(eta_fit, r3_data, B, N),
    r4 = ap_average_ml_prediction(eta_fit, r4_data, B, N)
  )
}

ap_canonical_validation_cache <- function(validation, fitted, learner,
                                          config, K, rep_id, fold_id) {
  n <- nrow(validation)
  ids <- sort(unique(validation$id))
  eta <- list(correct = list(), wrong = list())
  for (status in c("correct", "wrong")) {
    for (arm in 0:1) {
      eta[[status]][[as.character(arm)]] <- ap_predict_eta(
        fitted$eta[[status]][[as.character(arm)]], validation
      )
    }
  }

  theta_ratio <- list()
  tau_ratio <- list()
  clipping <- list()
  for (ks in c("correct", "wrong")) {
    theta_ratio[[ks]] <- list()
    tau_ratio[[ks]] <- list()
    clipping[[ks]] <- list()
    for (xs in c("correct", "wrong")) {
      theta_ratio[[ks]][[xs]] <- list(
        theta11 = numeric(n), theta10 = numeric(n), theta00 = numeric(n)
      )
      tau_ratio[[ks]][[xs]] <- numeric(n)
      clipping[[ks]][[xs]] <- list(
        theta11 = logical(length(ids)), theta10 = logical(length(ids)),
        theta00 = logical(length(ids)), tau = numeric(0)
      )
    }
  }

  theta_parts <- list()
  tau_parts <- list()
  for (es in c("correct", "wrong")) {
    theta_parts[[es]] <- list()
    tau_parts[[es]] <- list()
    for (ps in c("correct", "wrong")) {
      theta_parts[[es]][[ps]] <- list()
      tau_parts[[es]][[ps]] <- list()
      for (ks in c("correct", "wrong")) {
        theta_parts[[es]][[ps]][[ks]] <- lapply(
          c("theta11", "theta10", "theta00"),
          function(x) list(u1 = numeric(n), u2 = numeric(n), u3 = numeric(n))
        )
        names(theta_parts[[es]][[ps]][[ks]]) <-
          c("theta11", "theta10", "theta00")
        tau_parts[[es]][[ps]][[ks]] <- list(
          r1 = numeric(n), r2 = numeric(n), r3 = numeric(n), r4 = numeric(n)
        )
      }
    }
  }

  definitions <- list(
    theta11 = c(a = 1, a_prime = 1),
    theta10 = c(a = 1, a_prime = 0),
    theta00 = c(a = 0, a_prime = 0)
  )

  # Parametric partial means are analytic and vectorize over the whole fold.
  if (learner == "par") {
    for (es in c("correct", "wrong")) {
      for (ps in c("correct", "wrong")) {
        for (ks in c("correct", "wrong")) {
          for (component in names(definitions)) {
            definition <- definitions[[component]]
            theta_parts[[es]][[ps]][[ks]][[component]] <-
              ap_partial_theta_par(
                fitted$eta[[es]][[as.character(definition[["a"]])]],
                validation, fitted$phi[[ps]], fitted$kappa[[ks]],
                definition[["a"]], definition[["a_prime"]]
              )
          }
          tau_parts[[es]][[ps]][[ks]] <- ap_partial_tau_par(
            fitted$eta[[es]][["1"]], validation,
            fitted$phi[[ps]], fitted$kappa[[ks]]
          )
        }
      }
    }
  }

  for (id_index in seq_along(ids)) {
    id <- ids[[id_index]]
    index <- which(validation$id == id)
    cluster <- validation[index, , drop = FALSE]
    for (ks in c("correct", "wrong")) {
      for (xs in c("correct", "wrong")) {
        for (component in names(definitions)) {
          definition <- definitions[[component]]
          ratio <- ap_canonical_theta_ratio_cluster(
            cluster, fitted$kappa[[ks]], fitted$xi[[xs]],
            definition[["a"]], definition[["a_prime"]], config
          )
          theta_ratio[[ks]][[xs]][[component]][index] <- ratio$value
          clipping[[ks]][[xs]][[component]][id_index] <- ratio$clipped
        }
        ratio_tau <- ap_canonical_tau_ratio_cluster(
          cluster, fitted$kappa[[ks]], fitted$xi[[xs]], config
        )
        tau_ratio[[ks]][[xs]][index] <- ratio_tau$value
        clipping[[ks]][[xs]]$tau <- c(
          clipping[[ks]][[xs]]$tau, ratio_tau$clipped
        )
      }
    }

    draws <- NULL
    if (learner != "par") {
      integration_draws <- ap_ml_integration_draws(config, K)
      draws <- list(phi = list(), kappa = list())
      for (status in c("correct", "wrong")) {
        draws$phi[[status]] <- list()
        draws$kappa[[status]] <- list()
        for (arm in 0:1) {
          draws$phi[[status]][[as.character(arm)]] <- ap_draw_process(
            fitted$phi[[status]], cluster, arm, integration_draws,
            ap_hash_seed(config$base_seed, "canonical-mc", learner, K,
                         rep_id, fold_id, id, "phi", status, arm)
          )
          draws$kappa[[status]][[as.character(arm)]] <- ap_draw_process(
            fitted$kappa[[status]], cluster, arm, integration_draws,
            ap_hash_seed(config$base_seed, "canonical-mc", learner, K,
                         rep_id, fold_id, id, "kappa", status, arm)
          )
        }
      }
    }

    if (learner != "par") {
      for (es in c("correct", "wrong")) {
        for (ps in c("correct", "wrong")) {
          for (ks in c("correct", "wrong")) {
            for (component in names(definitions)) {
              definition <- definitions[[component]]
              part <- ap_partial_theta_ml(
                fitted$eta[[es]][[as.character(definition[["a"]])]], cluster,
                draws$kappa[[ks]][[as.character(definition[["a_prime"]])]],
                draws$phi[[ps]][[as.character(definition[["a"]])]]
              )
              for (name in names(part)) {
                theta_parts[[es]][[ps]][[ks]][[component]][[name]][index] <-
                  part[[name]]
              }
            }
            part_tau <- ap_partial_tau_ml(
              fitted$eta[[es]][["1"]], cluster,
              draws$kappa[[ks]][["0"]], draws$kappa[[ks]][["1"]],
              draws$phi[[ps]][["1"]]
            )
            for (name in names(part_tau)) {
              tau_parts[[es]][[ps]][[ks]][[name]][index] <- part_tau[[name]]
            }
          }
        }
      }
    }
  }

  list(
    eta = eta, theta_ratio = theta_ratio, tau_ratio = tau_ratio,
    theta_parts = theta_parts, tau_parts = tau_parts, clipping = clipping
  )
}

ap_canonical_clipping_diagnostics <- function(clipping, learner, K, rep_id) {
  method <- paste0("canonical_eif_", learner, "_stabilized")
  rows <- list()
  index <- 1L
  for (fold_id in seq_along(clipping)) {
    fold <- clipping[[fold_id]]
    for (ks in names(fold)) {
      for (xs in names(fold[[ks]])) {
        for (component in names(fold[[ks]][[xs]])) {
          indicator <- as.logical(fold[[ks]][[xs]][[component]])
          rows[[index]] <- data.frame(
            study = "canonical", K = K, rep = rep_id, method = method,
            fold = fold_id, scenario = NA_character_, component = component,
            quantity = "density_ratio",
            detail = paste0("kappa_", ks, ";xi_", xs),
            n = length(indicator), truncated = sum(indicator),
            truncated_fraction = mean(indicator)
          )
          index <- index + 1L
        }
      }
    }
  }
  do.call(rbind, rows)
}

ap_stabilization_epsilon <- function(Y, eta, clever) {
  ok <- is.finite(Y) & is.finite(eta) & is.finite(clever) & clever >= 0
  if (sum(ok) < 2L || sum(clever[ok]) <= 0) return(0)
  sum(clever[ok] * (Y[ok] - eta[ok])) / sum(clever[ok])
}

ap_score_theta_fold <- function(validation, a, a_prime, eta, ratio, parts,
                                config, stabilize = TRUE) {
  pi_a <- if (a == 1) config$dgp$treatment_prob else
    1 - config$dgp$treatment_prob
  pi_ap <- if (a_prime == 1) config$dgp$treatment_prob else
    1 - config$dgp$treatment_prob
  target <- validation$A == a
  epsilon_C <- epsilon_I <- 0
  if (stabilize) {
    epsilon_C <- ap_stabilization_epsilon(
      validation$Y[target], eta[target],
      ratio[target] / pi_a / validation$N[target]
    )
    epsilon_I <- ap_stabilization_epsilon(
      validation$Y[target], eta[target], ratio[target] / pi_a
    )
  }
  ids <- sort(unique(validation$id))
  out <- ap_empty_component(ids)
  score <- function(epsilon) {
    I1 <- as.numeric(validation$A == a) / pi_a * ratio *
      (validation$Y - eta - epsilon)
    I2 <- as.numeric(validation$A == a) / pi_a * (parts$u1 - parts$u3)
    I3 <- as.numeric(validation$A == a_prime) / pi_ap * (parts$u2 - parts$u3)
    I4 <- parts$u3 + epsilon
    I1 + I2 + I3 + I4
  }
  sum_C <- drop(rowsum(score(epsilon_C), validation$id, reorder = TRUE))
  sum_I <- drop(rowsum(score(epsilon_I), validation$id, reorder = TRUE))
  names(sum_C) <- names(sum_I) <- rownames(rowsum(
    rep(1, nrow(validation)), validation$id, reorder = TRUE
  ))
  N_by_id <- setNames(ap_cluster_frame(validation)$N,
                      as.character(ap_cluster_frame(validation)$id))
  out$C[names(sum_C)] <- sum_C / N_by_id[names(sum_C)]
  out$I[names(sum_I)] <- sum_I
  attr(out, "epsilon_C") <- epsilon_C
  attr(out, "epsilon_I") <- epsilon_I
  out
}

ap_score_tau_fold <- function(validation, eta, ratio, parts, config,
                              stabilize = TRUE) {
  pi <- config$dgp$treatment_prob
  target <- validation$A == 1
  epsilon_C <- epsilon_I <- 0
  if (stabilize) {
    epsilon_C <- ap_stabilization_epsilon(
      validation$Y[target], eta[target],
      ratio[target] / pi / validation$N[target]
    )
    epsilon_I <- ap_stabilization_epsilon(
      validation$Y[target], eta[target], ratio[target] / pi
    )
  }
  ids <- sort(unique(validation$id))
  out <- ap_empty_component(ids)
  score <- function(epsilon) {
    I1 <- validation$A / pi * ratio * (validation$Y - eta - epsilon)
    I2 <- validation$A / pi * (parts$r1 - parts$r4)
    I3 <- validation$A / pi * (parts$r2 - parts$r4)
    I4 <- (1 - validation$A) / (1 - pi) * (parts$r3 - parts$r4)
    I5 <- parts$r4 + epsilon
    I1 + I2 + I3 + I4 + I5
  }
  sum_C <- drop(rowsum(score(epsilon_C), validation$id, reorder = TRUE))
  sum_I <- drop(rowsum(score(epsilon_I), validation$id, reorder = TRUE))
  names(sum_C) <- names(sum_I) <- rownames(rowsum(
    rep(1, nrow(validation)), validation$id, reorder = TRUE
  ))
  N_by_id <- setNames(ap_cluster_frame(validation)$N,
                      as.character(ap_cluster_frame(validation)$id))
  out$C[names(sum_C)] <- sum_C / N_by_id[names(sum_C)]
  out$I[names(sum_I)] <- sum_I
  attr(out, "epsilon_C") <- epsilon_C
  attr(out, "epsilon_I") <- epsilon_I
  out
}

ap_score_canonical_cache <- function(validation, cache, scenario, config,
                                     stabilize = TRUE) {
  es <- ap_bool_name(scenario$eta)
  ps <- ap_bool_name(scenario$phi)
  ks <- ap_bool_name(scenario$kappa)
  xs <- ap_bool_name(scenario$xi)
  definitions <- list(
    theta11 = c(a = 1, a_prime = 1),
    theta10 = c(a = 1, a_prime = 0),
    theta00 = c(a = 0, a_prime = 0)
  )
  components <- list()
  for (component in names(definitions)) {
    definition <- definitions[[component]]
    components[[component]] <- ap_score_theta_fold(
      validation, definition[["a"]], definition[["a_prime"]],
      cache$eta[[es]][[as.character(definition[["a"]])]],
      cache$theta_ratio[[ks]][[xs]][[component]],
      cache$theta_parts[[es]][[ps]][[ks]][[component]], config, stabilize
    )
  }
  components$tau <- ap_score_tau_fold(
    validation, cache$eta[[es]][["1"]], cache$tau_ratio[[ks]][[xs]],
    cache$tau_parts[[es]][[ps]][[ks]], config, stabilize
  )
  components
}

ap_merge_fold_component <- function(target, source) {
  for (name in names(target)) {
    target[[name]]$C[names(source[[name]]$C)] <- source[[name]]$C
    target[[name]]$I[names(source[[name]]$I)] <- source[[name]]$I
  }
  target
}

ap_run_canonical_eif <- function(data, K, rep_id, learner, config,
                                 stabilize = TRUE) {
  scenarios <- ap_canonical_scenarios()
  ids <- sort(unique(data$id))
  components <- lapply(scenarios, function(x) {
    list(
      theta11 = ap_empty_component(ids), theta10 = ap_empty_component(ids),
      theta00 = ap_empty_component(ids), tau = ap_empty_component(ids)
    )
  })
  fold <- ap_make_folds(
    data, config$folds,
    ap_hash_seed(config$base_seed, "fold", K, rep_id)
  )
  diagnostics <- list()
  for (fold_id in seq_len(config$folds)) {
    validation_ids <- as.numeric(names(fold)[fold == fold_id])
    validation <- data[data$id %in% validation_ids, , drop = FALSE]
    train <- data[!data$id %in% validation_ids, , drop = FALSE]
    fitted <- ap_fit_canonical_nuisances(
      train, learner, config, K, rep_id, fold_id
    )
    cache <- ap_canonical_validation_cache(
      validation, fitted, learner, config, K, rep_id, fold_id
    )
    for (scenario_id in names(scenarios)) {
      scored <- ap_score_canonical_cache(
        validation, cache, scenarios[[scenario_id]], config, stabilize
      )
      components[[scenario_id]] <- ap_merge_fold_component(
        components[[scenario_id]], scored
      )
    }
    diagnostics[[fold_id]] <- cache$clipping
  }
  N_by_id <- setNames(ap_cluster_frame(data)$N, ap_cluster_frame(data)$id)
  method <- paste0("canonical_eif_", learner,
                   if (stabilize) "_stabilized" else "")
  rows <- lapply(names(scenarios), function(scenario_id) {
    effect_rows <- ap_effect_rows(components[[scenario_id]], N_by_id, config)
    ap_finalize_rows(
      effect_rows, "canonical", K, rep_id, scenario_id, method,
      scenarios[[scenario_id]]
    )
  })
  list(rows = do.call(rbind, rows), components = components,
       clipping = diagnostics,
       diagnostics = ap_canonical_clipping_diagnostics(
         diagnostics, learner, K, rep_id
       ),
       fold = fold)
}

ap_score_canonical_prefit <- function(validation, fitted, learner, config,
                                      K = nrow(ap_cluster_frame(validation)),
                                      rep_id = 0L, stabilize = TRUE) {
  cache <- ap_canonical_validation_cache(
    validation, fitted, learner, config, K, rep_id, 0L
  )
  ids <- sort(unique(validation$id))
  N_by_id <- setNames(ap_cluster_frame(validation)$N, ap_cluster_frame(validation)$id)
  output <- list()
  for (scenario in ap_canonical_scenarios()) {
    components <- ap_score_canonical_cache(
      validation, cache, scenario, config, stabilize
    )
    output[[scenario$id]] <- list(
      effects = ap_effect_rows(components, N_by_id, config),
      components = components
    )
  }
  list(output = output, clipping = cache$clipping)
}

ap_gcomp_components_par <- function(data, eta, phi, kappa) {
  ids <- sort(unique(data$id))
  components <- list(
    theta11 = ap_empty_component(ids), theta10 = ap_empty_component(ids),
    theta00 = ap_empty_component(ids), tau = ap_empty_component(ids)
  )
  definitions <- list(
    theta11 = c(a = 1, a_prime = 1),
    theta10 = c(a = 1, a_prime = 0),
    theta00 = c(a = 0, a_prime = 0)
  )
  cf <- ap_cluster_frame(data)
  N_by_id <- setNames(cf$N, as.character(cf$id))
  aggregate_value <- function(value) {
    sums <- drop(rowsum(value, data$id, reorder = TRUE))
    names(sums) <- rownames(rowsum(rep(1, nrow(data)), data$id, reorder = TRUE))
    list(C = sums / N_by_id[names(sums)], I = sums)
  }
  for (component in names(definitions)) {
    definition <- definitions[[component]]
    value <- ap_predict_eta_moment(
      eta[[as.character(definition[["a"]])]], data,
      ap_moments_theta(data, phi, kappa,
        definition[["a"]], definition[["a_prime"]])
    )
    aggregate <- aggregate_value(value)
    components[[component]]$C[names(aggregate$C)] <- aggregate$C
    components[[component]]$I[names(aggregate$I)] <- aggregate$I
  }
  value_tau <- ap_predict_eta_moment(
    eta[["1"]], data, ap_moments_tau(data, phi, kappa)
  )
  aggregate <- aggregate_value(value_tau)
  components$tau$C[names(aggregate$C)] <- aggregate$C
  components$tau$I[names(aggregate$I)] <- aggregate$I
  components
}

ap_fit_gcomp <- function(data, eta_correct, phi_correct, kappa_correct,
                         config, seed) {
  eta <- list()
  for (arm in 0:1) {
    eta[[as.character(arm)]] <- ap_fit_eta(
      data, arm, eta_correct, "par", seed + arm
    )
  }
  phi <- ap_fit_process(
    data, "Z", phi_correct, FALSE, "par", seed + 11L
  )
  if (kappa_correct) {
    phi_for_kappa <- if (phi_correct) phi else ap_fit_process(
      data, "Z", TRUE, FALSE, "par", seed + 19L
    )
    xi_for_kappa <- ap_fit_process(
      data, "M", TRUE, TRUE, "par", seed + 21L
    )
    kappa <- ap_derive_compatible_kappa(
      phi_for_kappa, xi_for_kappa, data
    )
  } else {
    kappa <- ap_fit_process(data, "M", FALSE, FALSE, "par", seed + 23L)
  }
  list(eta = eta, phi = phi, kappa = kappa)
}

ap_gcomp_point <- function(data, scenario, config, seed) {
  fitted <- ap_fit_gcomp(
    data, scenario$eta, scenario$phi, scenario$kappa, config, seed
  )
  components <- ap_gcomp_components_par(
    data, fitted$eta, fitted$phi, fitted$kappa
  )
  N_by_id <- setNames(ap_cluster_frame(data)$N, ap_cluster_frame(data)$id)
  ap_effect_rows(components, N_by_id, config)[, c("estimand", "effect", "scale", "point")]
}

ap_gcomp_point_matrix <- function(data, config, seed) {
  eta <- list(correct = list(), wrong = list())
  for (correct in c(TRUE, FALSE)) {
    status <- ap_bool_name(correct)
    for (arm in 0:1) {
      eta[[status]][[as.character(arm)]] <- ap_fit_eta(
        data, arm, correct, "par", seed + 10L * as.integer(!correct) + arm
      )
    }
  }
  phi <- list(
    correct = ap_fit_process(data, "Z", TRUE, FALSE, "par", seed + 101L),
    wrong = ap_fit_process(data, "Z", FALSE, FALSE, "par", seed + 102L)
  )
  xi_for_kappa <- ap_fit_process(
    data, "M", TRUE, TRUE, "par", seed + 201L
  )
  kappa <- list(
    correct = ap_derive_compatible_kappa(
      phi$correct, xi_for_kappa, data
    ),
    wrong = ap_fit_process(data, "M", FALSE, FALSE, "par", seed + 202L)
  )
  combinations <- expand.grid(
    eta = c(TRUE, FALSE), phi = c(TRUE, FALSE), kappa = c(TRUE, FALSE),
    KEEP.OUT.ATTRS = FALSE
  )
  output <- matrix(
    NA_real_, nrow(combinations), length(ap_truth(config)),
    dimnames = list(NULL, names(ap_truth(config)))
  )
  combo_id <- character(nrow(combinations))
  cf <- ap_cluster_frame(data)
  N_by_id <- setNames(cf$N, cf$id)
  for (i in seq_len(nrow(combinations))) {
    combo <- combinations[i, ]
    es <- ap_bool_name(combo$eta)
    ps <- ap_bool_name(combo$phi)
    ks <- ap_bool_name(combo$kappa)
    components <- ap_gcomp_components_par(
      data, eta[[es]], phi[[ps]], kappa[[ks]]
    )
    rows <- ap_effect_rows(components, N_by_id, config)
    output[i, rows$estimand] <- rows$point
    combo_id[[i]] <- paste0(as.integer(c(combo$eta, combo$phi, combo$kappa)),
                            collapse = "")
  }
  rownames(output) <- combo_id
  output
}

ap_cluster_bootstrap <- function(data, seed) {
  ids <- sort(unique(data$id))
  set.seed(seed)
  sampled <- sample(ids, length(ids), replace = TRUE)
  pieces <- vector("list", length(sampled))
  for (i in seq_along(sampled)) {
    piece <- data[data$id == sampled[[i]], , drop = FALSE]
    piece$id <- i
    pieces[[i]] <- piece
  }
  do.call(rbind, pieces)
}

ap_run_gcomp <- function(data, K, rep_id, config) {
  scenarios <- ap_canonical_scenarios()
  point <- ap_gcomp_point_matrix(
    data, config, ap_hash_seed(config$base_seed, "gcomp", K, rep_id)
  )
  boot <- array(
    NA_real_, c(config$gcomp_boot, nrow(point), ncol(point)),
    dimnames = list(NULL, rownames(point), colnames(point))
  )
  for (b in seq_len(config$gcomp_boot)) {
    boot_data <- ap_cluster_bootstrap(
      data, ap_hash_seed(config$base_seed, "gcomp-bootstrap", K, rep_id, b)
    )
    boot[b, , ] <- ap_gcomp_point_matrix(
      boot_data, config,
      ap_hash_seed(config$base_seed, "gcomp-fit", K, rep_id, b)
    )
  }
  rows <- list()
  row_index <- 1L
  truth <- ap_truth(config)
  for (scenario in scenarios) {
    combo <- paste0(as.integer(c(scenario$eta, scenario$phi, scenario$kappa)),
                    collapse = "")
    for (target in names(truth)) {
      estimate <- point[combo, target]
      boot_values <- boot[, combo, target]
      interval <- quantile(
        boot_values, c(0.025, 0.975), na.rm = TRUE, names = FALSE
      )
      effect <- sub("_[CI]$", "", target)
      scale <- if (endsWith(target, "_C")) "cluster" else "individual"
      one <- data.frame(
        estimand = target, effect = effect, scale = scale, point = estimate,
        se = sd(boot_values, na.rm = TRUE), ci_low = interval[[1L]],
        ci_up = interval[[2L]], true = unname(truth[[target]])
      )
      theory <- scenario
      theory$point_consistency_expected <- scenario$gcomp_point_expected
      theory$asymptotic_linearity_expected <- scenario$gcomp_point_expected
      theory$ci_validity_expected <- scenario$gcomp_point_expected
      rows[[row_index]] <- ap_finalize_rows(
        one, "canonical", K, rep_id, scenario$id, "gcomp_parametric",
        theory, ci_method = "cluster_percentile_bootstrap"
      )
      row_index <- row_index + 1L
    }
  }
  do.call(rbind, rows)
}
