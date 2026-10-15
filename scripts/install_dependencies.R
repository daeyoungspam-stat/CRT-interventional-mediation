options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 1L)
packages <- c("SuperLearner", "glmnet", "ranger", "data.table", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, dependencies = c("Depends", "Imports", "LinkingTo"))
}
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
print(data.frame(package = packages, version = vapply(
  packages, function(x) as.character(packageVersion(x)), ""
)))

