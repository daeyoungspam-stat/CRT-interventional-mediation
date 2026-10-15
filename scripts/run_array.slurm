#!/usr/bin/env bash
#SBATCH --job-name=crtmedsim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
# Set R_MODULE to the module provided by the execution environment.
if [[ -n "${R_MODULE:-}" ]]; then module load "$R_MODULE"; fi
Rscript scripts/array.R "$@"

