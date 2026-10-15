# Slurm execution

Install the package and dependencies in the R environment available on compute
nodes. When modules are available, select the appropriate module in the allocated
environment. The array script accepts `R_MODULE` and loads it on each compute node.

For a project-specific package library, run in that R environment:

```sh
export R_LIBS_USER="$HOME/Rlibs/crtmedsim"
mkdir -p "$R_LIBS_USER"
Rscript scripts/install_dependencies.R
R CMD INSTALL -l "$R_LIBS_USER" .
```

Retain `R_LIBS_USER` in the submission environment so every task uses that library.

From the repository root, create the workspace once:

```sh
mkdir -p logs
Rscript -e 'library(crtmedsim); create_simulation("runs/full", simulation_config("full"))'
```

Submit both studies. Each array has 800 tasks: 200 chunks of 10 replications for
each of four cluster counts. The following commands allow 20 simultaneous tasks
per study with 10 replication workers per task; scheduler limits determine actual
concurrency. Supply partition/account options required by the execution environment.

```sh
sbatch --array=1-800%20 scripts/run_array.slurm runs/full canonical 10
sbatch --array=1-800%20 scripts/run_array.slurm runs/full lowdim 10
```

Keep the repository as the submission working directory. The workspace must be
visible on every allocated node. Resubmit an interrupted array or selected task
IDs with the same arguments to reuse completed matching replications.

After both arrays finish, run aggregation in an allocated R environment:

```sh
Rscript -e 'library(crtmedsim); collect_results("runs/full")'
```

The combined bundle is `runs/full/results/annals_simulation_bundle.rds`.
The `logs/` directory contains scheduler logs; workspace `logs/` contains engine
logs. Adjust time, CPUs, concurrency, and partition options to the available
allocation. Changing the chunk size requires recalculating the array length as
`length(K_values) * ceiling(reps / chunk_size)`.

