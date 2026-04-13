#!/usr/bin/env bash
#SBATCH --job-name=axi-ex03-3d
#SBATCH --partition=windfall
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --chdir=/groups/wolf_lab/jon/AxiSEM3D/examples/03_salt_body_SEG_local
#SBATCH --output=/groups/wolf_lab/jon/AxiSEM3D/examples/03_salt_body_SEG_local/slurm-%j.out
#SBATCH --error=/groups/wolf_lab/jon/AxiSEM3D/examples/03_salt_body_SEG_local/slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/simu3D/output_windfall_${SLURM_JOB_ID}"

./run3D.sh