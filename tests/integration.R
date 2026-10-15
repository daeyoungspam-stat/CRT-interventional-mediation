if (identical(Sys.getenv("CRTMEDSIM_INTEGRATION"), "true")) {
  library(crtmedsim)
  workspace <- tempfile("integration-")
  create_simulation(workspace, simulation_config("smoke"))
  run_simulation(workspace, cores = 1L)
  bundle <- collect_results(workspace)
  result <- read_results(bundle)
  stopifnot(nrow(result$raw) == 684L, nrow(result$summary) == 684L,
            all(result$summary$n_rep == 1L),
            all(is.finite(result$raw$point)), all(is.finite(result$raw$se)),
            all(result$raw$ci_low <= result$raw$ci_up))
  files <- list.files(file.path(workspace, "results", "parts"),
                       recursive = TRUE, full.names = TRUE)
  before <- tools::md5sum(files)
  run_simulation(workspace, cores = 1L)
  stopifnot(identical(before, tools::md5sum(files)))
  report <- file.path(workspace, "results", "report")
  stopifnot(file.exists(file.path(report, "simulation_summary.csv")),
            file.exists(file.path(report, "comprehensive_simulation_tables.tex")),
            length(list.files(file.path(report, "figures"), pattern = "\\.pdf$")) > 0L)
  unlink(workspace, recursive = TRUE)
}
