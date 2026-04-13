#!/usr/bin/env bash
#SBATCH --job-name=axi-ex04-plume
#SBATCH --partition=windfall
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --chdir=/groups/wolf_lab/jon/AxiSEM3D/examples/04_simple_3d_shapes/example_single_plume
#SBATCH --output=/groups/wolf_lab/jon/AxiSEM3D/examples/04_simple_3d_shapes/example_single_plume/slurm-%j.out
#SBATCH --error=/groups/wolf_lab/jon/AxiSEM3D/examples/04_simple_3d_shapes/example_single_plume/slurm-%j.err

set -euo pipefail

export NP="${NP:-4}"
export OUTPUT_DIR="$PWD/output_windfall_${SLURM_JOB_ID}"

./run.sh