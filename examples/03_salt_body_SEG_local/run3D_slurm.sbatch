#!/usr/bin/env bash
#SBATCH --job-name=axi-ex03-3d
#SBATCH --partition=FIXME
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/simu3D/output_slurm_${SLURM_JOB_ID}"

./run3D.sh
