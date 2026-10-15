#!/usr/bin/env Rscript

source("annals_core.R")

ap_assert(requireNamespace("data.table", quietly = TRUE),
          "data.table is required.")
ap_assert(requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 is required.")

args <- ap_parse_args(commandArgs(trailingOnly = TRUE))
bundle_path <- ap_arg(
  args, "bundle", file.path("final_results", "annals_simulation_bundle.rds")
)
output_dir <- ap_arg(args, "output-dir", file.path("final_results", "report"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(bundle_path)
ap_assert(isTRUE(bundle$complete), "The result bundle is incomplete.")
ap_assert(identical(bundle$code_signature, ap_estimator_code_signature()),
          "The result bundle and estimator source code differ.")
require_certification <- identical(
  tolower(ap_arg(args, "require-certification", "true")), "true"
)
if (require_certification) {
  ap_assert(!is.null(bundle$design_certification) &&
    isTRUE(bundle$design_certification$passed),
    "The result bundle lacks a passed design certification.")
}
raw <- data.table::as.data.table(bundle$raw)
summary <- data.table::as.data.table(bundle$summary)

canonical_order <- names(ap_canonical_scenarios())
lowdim_order <- names(ap_lowdim_scenarios())
method_order <- c(
  "canonical_eif_par_stabilized", "canonical_eif_ml_stabilized",
  "gcomp_parametric",
  "ld_gaussian_glmnet_stabilized", "ld_gaussian_ranger_stabilized",
  "ld_gaussian_sl_stabilized", "ld_resid_boot_glmnet_stabilized",
  "ld_resid_boot_ranger_stabilized", "ld_resid_boot_sl_stabilized"
)
palette <- c(
  "Canonical EIF: parametric" = "#287271",
  "Canonical EIF: machine learning" = "#E07A5F",
  "Parametric g-computation" = "#3D405B",
  "LD Gaussian + glmnet" = "#2A6FBB",
  "LD Gaussian + random forest" = "#7A5195",
  "LD Gaussian + Super Learner" = "#2F8F5B",
  "LD residual bootstrap + glmnet" = "#D1495B",
  "LD residual bootstrap + random forest" = "#B07D2B",
  "LD residual bootstrap + Super Learner" = "#5C677D"
)

base_theme <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    axis.text.x = ggplot2::element_text(angle = 50, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

ap_plot_pair <- function(study, K, effect, scale, scenarios = NULL,
                         methods = NULL, prefix = "all") {
  study_value <- study
  K_value <- K
  effect_value <- effect
  scale_value <- scale
  data <- raw[
    study == study_value & K == K_value & effect == effect_value &
      scale == scale_value
  ]
  sums <- summary[
    study == study_value & K == K_value & effect == effect_value &
      scale == scale_value
  ]
  if (!is.null(scenarios)) {
    data <- data[scenario %in% scenarios]
    sums <- sums[scenario %in% scenarios]
  }
  if (!is.null(methods)) {
    data <- data[method %in% methods]
    sums <- sums[method %in% methods]
  }
  if (!nrow(data)) return(character(0))
  scenario_levels <- if (study == "canonical") canonical_order else lowdim_order
  scenario_levels <- scenario_levels[scenario_levels %in% unique(data$scenario)]
  data[, scenario := factor(scenario, levels = scenario_levels)]
  sums[, scenario := factor(scenario, levels = scenario_levels)]
  data[, method := factor(method, levels = method_order)]
  sums[, method := factor(method, levels = method_order)]
  data[, method_label := factor(method_label, levels = names(palette))]
  sums[, method_label := factor(method_label, levels = names(palette))]
  truth <- unique(data$true)
  truth <- truth[[1L]]
  study_label <- if (study == "canonical") "Canonical" else "Low-dimensional"

  estimation <- ggplot2::ggplot(
    data,
    ggplot2::aes(scenario, point, fill = method_label)
  ) +
    ggplot2::geom_boxplot(
      width = 0.72, outlier.shape = NA, position = ggplot2::position_dodge(0.78)
    ) +
    ggplot2::geom_hline(
      yintercept = truth, color = "#B22222", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = TRUE) +
    ggplot2::labs(
      x = NULL, y = "Estimate", fill = NULL,
      title = sprintf("%s, K=%d, %s (%s scale)",
                      study_label, K, effect, scale)
    ) + base_theme

  coverage_lower <- max(0, min(0.85, min(sums$coverage, na.rm = TRUE) - 0.03))
  coverage <- ggplot2::ggplot(
    sums,
    ggplot2::aes(scenario, coverage, color = method_label, group = method_label)
  ) +
    ggplot2::geom_hline(
      yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_line(
      position = ggplot2::position_dodge(width = 0.30), linewidth = 0.45
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.30), size = 1.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, coverage - 1.96 * coverage_mcse),
        ymax = pmin(1, coverage + 1.96 * coverage_mcse)
      ),
      position = ggplot2::position_dodge(width = 0.30), width = 0.12,
      linewidth = 0.35
    ) +
    ggplot2::coord_cartesian(ylim = c(coverage_lower, 1.002)) +
    ggplot2::scale_color_manual(values = palette, drop = TRUE) +
    ggplot2::labs(
      x = NULL, y = "Coverage", color = NULL,
      title = sprintf("Coverage: %s, K=%d, %s (%s scale)",
                      study_label, K, effect, scale)
    ) + base_theme

  stem <- sprintf("%s_%s_K%d_%s_%s", prefix, study, K, effect, scale)
  estimation_path <- file.path(figure_dir, paste0(stem, "_estimation.pdf"))
  coverage_path <- file.path(figure_dir, paste0(stem, "_coverage.pdf"))
  ggplot2::ggsave(estimation_path, estimation, width = 9.0, height = 4.5,
                  device = grDevices::pdf)
  ggplot2::ggsave(coverage_path, coverage, width = 9.0, height = 4.5,
                  device = grDevices::pdf)
  created <- c(estimation_path, coverage_path)
  ap_assert(all(file.exists(created)) && all(file.info(created)$size > 0),
            paste("A PDF figure was not created for", stem))
  created
}

all_figures <- character(0)
for (study in c("canonical", "lowdim")) {
  for (K in bundle$config$K_values) {
    for (effect in c("IDE", "IIME", "ISME")) {
      for (scale in c("cluster", "individual")) {
        all_figures <- c(all_figures, ap_plot_pair(
          study, K, effect, scale, prefix = "appendix"
        ))
      }
    }
  }
}

canonical_main_scenarios <- c(
  "C-1111", "C-0111", "C-1011", "C-1101", "C-1110", "C-0000"
)
canonical_main_methods <- c(
  "canonical_eif_par_stabilized", "canonical_eif_ml_stabilized",
  "gcomp_parametric"
)
lowdim_main_scenarios <- c("LD-AC", paste0("LD-T", 2:6))
lowdim_main_methods <- c(
  "ld_gaussian_sl_stabilized", "ld_resid_boot_sl_stabilized"
)
main_figures <- character(0)
for (K in bundle$config$K_values) {
  for (effect in c("IDE", "ISME")) {
    main_figures <- c(main_figures, ap_plot_pair(
      "canonical", K, effect, "cluster", canonical_main_scenarios,
      canonical_main_methods, "main"
    ))
    main_figures <- c(main_figures, ap_plot_pair(
      "lowdim", K, effect, "cluster", lowdim_main_scenarios,
      lowdim_main_methods, "main"
    ))
  }
}

ap_latex_escape <- function(x) {
  x <- gsub("_", "\\_", as.character(x), fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x
}

scenario_file <- file.path(output_dir, "scenario_guarantees.tex")
connection <- file(scenario_file, open = "wt")
writeLines(c(
  "% Auto-generated scenario guarantees.",
  "\\noindent\\textit{Note.} P, AL, and CI denote point consistency, asymptotic linearity, and interval validity under the corresponding untruncated theoretical nuisance limits. The implementation uses only the prespecified fixed probability bounds $[0.025,0.975]$ and ratio bounds $[0.05,20]$.",
  "\\begingroup\\scriptsize",
  "\\begin{longtable}{l*{13}{c}}",
  "\\caption{Canonical nuisance scenarios: point consistency, asymptotic linearity, and confidence-interval validity.}\\\\",
  "\\toprule",
  "Scenario & $\\eta$ & $\\phi$ & $\\kappa$ & $\\xi$ & IDE-P & IIME-P & ISME-P & IDE-AL & IIME-AL & ISME-AL & IDE-CI & IIME-CI & ISME-CI \\\\",
  "\\midrule"
), connection)
for (i in seq_len(nrow(bundle$canonical_scenarios))) {
  x <- bundle$canonical_scenarios[i]
  writeLines(sprintf(
    "%s & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\",
    ap_latex_escape(x$scenario), x$eta, x$phi, x$kappa, x$xi,
    x$IDE_point, x$IIME_point, x$ISME_point,
    x$IDE_AL, x$IIME_AL, x$ISME_AL,
    x$IDE_CI, x$IIME_CI, x$ISME_CI
  ), connection)
}
writeLines(c("\\bottomrule", "\\end{longtable}", "\\endgroup"), connection)
writeLines(c(
  "\\begingroup\\scriptsize",
  "\\begin{longtable}{l*{13}{c}}",
  "\\caption{Parametric g-computation scenarios: point consistency, asymptotic linearity, and percentile-interval validity.}\\\\",
  "\\toprule",
  "Scenario & $\\eta$ & $\\phi$ & $\\kappa$ & $\\xi$ & IDE-P & IIME-P & ISME-P & IDE-AL & IIME-AL & ISME-AL & IDE-CI & IIME-CI & ISME-CI \\\\",
  "\\midrule"
), connection)
for (i in seq_len(nrow(bundle$gcomp_scenarios))) {
  x <- bundle$gcomp_scenarios[i]
  writeLines(sprintf(
    "%s & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\",
    ap_latex_escape(x$scenario), x$eta, x$phi, x$kappa, x$xi,
    x$IDE_point, x$IIME_point, x$ISME_point,
    x$IDE_AL, x$IIME_AL, x$ISME_AL,
    x$IDE_CI, x$IIME_CI, x$ISME_CI
  ), connection)
}
writeLines(c("\\bottomrule", "\\end{longtable}", "\\endgroup"), connection)
writeLines(c(
  "\\begingroup\\scriptsize",
  "\\begin{longtable}{l*{15}{c}}",
  "\\caption{Low-dimensional scenarios: nuisance status and separate inferential guarantees.}\\\\",
  "\\toprule",
  "Scenario & $\\eta_\\theta$ & $d$ & $s$ & $\\eta_\\tau$ & $e$ & $s_{-j}$ & IDE-P & IIME-P & ISME-P & IDE-AL & IIME-AL & ISME-AL & IDE-CI & IIME-CI & ISME-CI \\\\",
  "\\midrule"
), connection)
for (i in seq_len(nrow(bundle$lowdim_scenarios))) {
  x <- bundle$lowdim_scenarios[i]
  writeLines(sprintf(
    "%s & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\",
    ap_latex_escape(x$scenario), x$theta_eta, x$d, x$s,
    x$tau_eta, x$e, x$s_minus,
    x$IDE_point, x$IIME_point, x$ISME_point,
    x$IDE_AL, x$IIME_AL, x$ISME_AL,
    x$IDE_CI, x$IIME_CI, x$ISME_CI
  ), connection)
}
writeLines(c("\\bottomrule", "\\end{longtable}", "\\endgroup"), connection)
close(connection)

table_file <- file.path(output_dir, "comprehensive_simulation_tables.tex")
connection <- file(table_file, open = "wt")
writeLines(c(
  "% Auto-generated comprehensive numerical results.",
  "\\begingroup\\scriptsize",
  "\\begin{longtable}{lllrrlrrrrrr}",
  "\\caption{Monte Carlo performance across all studies, scenarios, methods, targets, and cluster counts.}\\\\",
  "\\toprule",
  "Study & Scenario & Method & $K$ & Target & Scale & Bias & ESD & ASE & RMSE & Coverage & Width \\\\",
  "\\midrule"
), connection)
ordered_summary <- summary[order(study, K, scenario, method, estimand)]
for (i in seq_len(nrow(ordered_summary))) {
  x <- ordered_summary[i]
  writeLines(sprintf(
    "%s & %s & %s & %d & %s & %s & %.4f & %.4f & %.4f & %.4f & %.3f & %.4f \\\\",
    ap_latex_escape(x$study), ap_latex_escape(x$scenario),
    ap_latex_escape(x$method_label), x$K, ap_latex_escape(x$effect),
    ap_latex_escape(x$scale), x$bias, x$empirical_sd, x$mean_se,
    x$rmse, x$coverage, x$mean_ci_width
  ), connection)
}
writeLines(c("\\bottomrule", "\\end{longtable}", "\\endgroup"), connection)
close(connection)

ap_write_figure_includes <- function(paths, path) {
  connection <- file(path, open = "wt")
  writeLines("% Auto-generated figure includes.", connection)
  for (figure in paths) {
    relative <- file.path("figures", basename(figure))
    writeLines(c(
      "\\begin{figure}[p]",
      "\\centering",
      sprintf("\\includegraphics[width=\\textwidth]{%s}", relative),
      sprintf("\\caption{%s}", ap_latex_escape(tools::file_path_sans_ext(
        basename(figure)
      ))),
      "\\end{figure}"
    ), connection)
  }
  close(connection)
}

ap_write_figure_includes(
  main_figures, file.path(output_dir, "main_figure_includes.tex")
)
ap_write_figure_includes(
  all_figures, file.path(output_dir, "appendix_figure_includes.tex")
)

cat("Report regenerated from one bundle: ", bundle_path, "\n", sep = "")
