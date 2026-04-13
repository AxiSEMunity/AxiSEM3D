#!/usr/bin/env bash
#SBATCH --job-name=axi-ex04-paper
#SBATCH --partition=windfall
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --chdir=/groups/wolf_lab/jon/AxiSEM3D/examples/04_simple_3d_shapes/example_release_paper
#SBATCH --output=/groups/wolf_lab/jon/AxiSEM3D/examples/04_simple_3d_shapes/example_release_paper/slurm-%j.out
#SBATCH --error=/groups/wolf_lab/jon/AxiSEM3D/examples/04_simple_3d_shapes/example_release_paper/slurm-%j.err

set -euo pipefail

export NP="${NP:-16}"
export OUTPUT_DIR="$PWD/output_windfall_${SLURM_JOB_ID}"

./run.sh