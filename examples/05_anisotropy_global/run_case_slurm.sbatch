#!/usr/bin/env bash
#SBATCH --partition=FIXME
#SBATCH --job-name=axi-ex05
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CASE_REL="${CASE_REL:?Set CASE_REL to a subcase path under examples/05_anisotropy_global}"
CASE_DIR="$SCRIPT_DIR/$CASE_REL"

if [[ ! -d "$CASE_DIR" ]]; then
  echo "Case directory not found: $CASE_REL" >&2
  exit 1
fi

export NP="${NP:-${SLURM_NTASKS:-16}}"
export OUTPUT_DIR="${OUTPUT_DIR:-$CASE_DIR/output_slurm_${SLURM_JOB_ID}}"

cd "$SCRIPT_DIR"
bash "$SCRIPT_DIR/run_case.sh" "$CASE_REL"
