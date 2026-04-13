#!/usr/bin/env bash
#SBATCH --partition=windfall
#SBATCH --job-name=axi-ex02-3d
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
RUN_DIR="$SCRIPT_DIR/simu_with_3d_crust"

export NP="${NP:-${SLURM_NTASKS:-4}}"
export OUTPUT_DIR="${OUTPUT_DIR:-$RUN_DIR/output_windfall_${SLURM_JOB_ID}}"

cd "$SCRIPT_DIR"
bash "$SCRIPT_DIR/run_with_3d_crust.sh"