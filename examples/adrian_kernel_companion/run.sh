#!/bin/bash
# Adrian Kernel Companion – run script
#
# Usage:
#   bash run.sh forward    – set up and run the forward simulation
#   bash run.sh backward   – set up and run the adjoint simulation
#
# Run forward first, then backward.  Then compute kernels with:
#   python compute_kernels.py
#
# Prerequisites for this script:
#   - axisem3d binary must be in this directory or on PATH
#   - MPI launcher available as mpirun
#
# axikernels itself is only needed later for kernel post-processing.

set -e

MODE="${1:-}"

if [[ "$MODE" != "forward" && "$MODE" != "backward" ]]; then
    echo "Usage: bash run.sh forward|backward"
    exit 1
fi

SIMU_DIR="simu_${MODE}"
INPUT_DIR="input_${MODE}"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory '$INPUT_DIR' not found."
    echo "       Input files are populated in Phase 2 of the setup plan."
    exit 1
fi

if [[ -f "$INPUT_DIR/PLACEHOLDER.txt" ]]; then
    echo "ERROR: $INPUT_DIR still contains placeholder content only."
    echo "       Populate the real AxiSEM3D input files in Phase 2 before running."
    exit 1
fi

if ! compgen -G "$INPUT_DIR/inparam*.yaml" > /dev/null; then
    echo "ERROR: no inparam YAML files were found in $INPUT_DIR."
    echo "       Expected AxiSEM3D input files are not present yet."
    exit 1
fi

# create simulation dir
mkdir -p "${SIMU_DIR}/input"

# copy input files
cp -r "${INPUT_DIR}/"* "${SIMU_DIR}/input/"

# locate binary (directory-local copy preferred, then PATH)
if [[ -f "./axisem3d" ]]; then
    cp ./axisem3d "${SIMU_DIR}/"
elif command -v axisem3d &>/dev/null; then
    cp "$(command -v axisem3d)" "${SIMU_DIR}/"
else
    echo "ERROR: axisem3d binary not found in this directory or on PATH."
    exit 1
fi

# run
cd "${SIMU_DIR}"
mpirun -np 4 ./axisem3d
cd ..

echo "Done. Simulation output is in ${SIMU_DIR}/output/"
