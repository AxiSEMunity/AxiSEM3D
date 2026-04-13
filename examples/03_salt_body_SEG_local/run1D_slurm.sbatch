#!/usr/bin/env bash
#SBATCH --job-name=axi-ex03-1d
#SBATCH --partition=FIXME
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/simu1D/output_slurm_${SLURM_JOB_ID}"

./run1D.sh
