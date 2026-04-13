#!/usr/bin/env bash
#SBATCH --job-name=axi-ex02-1d
#SBATCH --partition=windfall
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --chdir=/groups/wolf_lab/jon/AxiSEM3D/examples/02_3d_crust_S362ANI_regional
#SBATCH --output=/groups/wolf_lab/jon/AxiSEM3D/examples/02_3d_crust_S362ANI_regional/slurm-%j.out
#SBATCH --error=/groups/wolf_lab/jon/AxiSEM3D/examples/02_3d_crust_S362ANI_regional/slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/simu_with_1d_crust/output_windfall_${SLURM_JOB_ID}"

./run_with_1d_crust.sh