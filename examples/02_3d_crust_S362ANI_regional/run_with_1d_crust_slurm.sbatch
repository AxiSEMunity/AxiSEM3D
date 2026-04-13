#!/usr/bin/env bash
#SBATCH --job-name=axi-ex02-1d
#SBATCH --partition=FIXME
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/simu_with_1d_crust/output_slurm_${SLURM_JOB_ID}"

./run_with_1d_crust.sh
