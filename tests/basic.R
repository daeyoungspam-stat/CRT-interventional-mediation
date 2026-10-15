library(crtmedsim)

full <- simulation_config()
smoke <- simulation_config("smoke")
stopifnot(full$reps == 2000L, identical(full$K_values, c(50L, 75L, 100L, 150L)),
          full$folds == 5L, smoke$reps == 1L, identical(smoke$K_values, 50L))
for (name in setdiff(names(full), c("reps", "K_values"))) {
  stopifnot(identical(full[[name]], smoke[[name]]))
}
stopifnot(identical(full$canonical_ml_library, c("SL.glm", "SL.ranger")))
engine <- system.file("engine", package = "crtmedsim", mustWork = TRUE)
manifest <- read.csv(file.path(engine, "checksums.csv"), stringsAsFactors = FALSE)
stopifnot(identical(unname(tools::md5sum(file.path(engine, manifest$file))), manifest$md5))

set.seed(123)
seed <- .Random.seed
first <- simulate_trial(50L, 1L)
second <- simulate_trial(50L, 1L)
third <- simulate_trial(50L, 2L)
stopifnot(identical(first, second), !identical(first, third), identical(seed, .Random.seed),
          length(unique(first$id)) == 50L, all(first$N >= 10L & first$N <= 20L))
RNGkind("L'Ecuyer-CMRG")
kind <- RNGkind()
seed <- .Random.seed
stopifnot(identical(first, simulate_trial(50L, 1L)),
          identical(kind, RNGkind()), identical(seed, .Random.seed))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
expected <- c(IDE_C = 0.468017950368579, IIME_C = 0.330550604999998,
              ISME_C = 0.192445875000001, IDE_I = 0.489481156033335,
              IIME_I = 0.342964557, ISME_I = 0.188840114999999)
stopifnot(max(abs(simulation_truth() - expected)) < 1e-12)

path <- tempfile("workspace-")
create_simulation(path, smoke)
create_simulation(path, smoke)
stopifnot(identical(readRDS(file.path(path, "simulation_config.rds")), smoke))
changed <- smoke
changed$base_seed <- changed$base_seed + 1L
stopifnot(inherits(try(create_simulation(path, changed), silent = TRUE), "try-error"))
stopifnot(inherits(try(run_simulation(path, cores = 0), silent = TRUE), "try-error"))
stopifnot(inherits(try(run_simulation(path, rep_start = 2, rep_end = 1), silent = TRUE),
                   "try-error"))
unlink(path, recursive = TRUE)

occupied <- tempfile("occupied-")
dir.create(occupied)
writeLines("existing", file.path(occupied, "existing.txt"))
stopifnot(inherits(try(create_simulation(occupied, smoke), silent = TRUE), "try-error"))
stopifnot(identical(readLines(file.path(occupied, "existing.txt")), "existing"))
unlink(occupied, recursive = TRUE)

