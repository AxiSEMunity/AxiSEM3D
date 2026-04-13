#!/usr/bin/env bash
#SBATCH --job-name=axi-ex04-paper
#SBATCH --partition=FIXME
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

export NP="${NP:-16}"
export OUTPUT_DIR="$PWD/output_slurm_${SLURM_JOB_ID}"

./run.sh
