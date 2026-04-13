#!/usr/bin/env bash
#SBATCH --job-name=axi-ex01
#SBATCH --partition=FIXME
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/output_slurm_${SLURM_JOB_ID}"

./run.sh
