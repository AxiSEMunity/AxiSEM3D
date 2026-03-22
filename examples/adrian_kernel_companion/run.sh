#!/bin/bash
# Adrian Kernel Companion – run script
#
# Usage:
#   bash run.sh [forward]          – set up and run the forward simulation
#
#   The backward (adjoint) simulation is generated and run automatically by
#   axikernels when you call:
#       python compute_kernels.py
#   You do NOT need to run "bash run.sh backward" manually.
#
# Environment variables:
#   NRANKS         – number of MPI ranks to use (default: 4)
#   CLEAN_FORWARD  – set to 1 to remove an existing simu_forward/ first
#
# Prerequisites:
#   - axisem3d binary must be in this directory or on PATH
#   - axikernels installed: pip install -e /path/to/AxiSEM3D_Kernels
#
# The MPI launcher is auto-detected from the linked libraries of the axisem3d
# binary (mirrors the same logic used by axikernels for the backward run).

set -euo pipefail

# ── configuration ────────────────────────────────────────────────────────────
NRANKS="${NRANKS:-8}"
MODE="${1:-forward}"

if [[ "$MODE" != "forward" ]]; then
    echo "Usage: bash run.sh [forward]"
    echo ""
    echo "The backward simulation is generated automatically by axikernels."
    echo "Run 'python compute_kernels.py' after the forward simulation completes."
    exit 1
fi

SIMU_DIR="simu_forward"
INPUT_DIR="input_forward"

if [[ -d "$SIMU_DIR/output" || -f "$SIMU_DIR/axisem3d" ]]; then
    if [[ "${CLEAN_FORWARD:-0}" == "1" ]]; then
        echo "Removing existing ${SIMU_DIR}/ because CLEAN_FORWARD=1"
        rm -rf "$SIMU_DIR"
    else
        echo "ERROR: ${SIMU_DIR}/ already appears to contain a previous run."
        echo "       Remove it manually, or rerun with: CLEAN_FORWARD=1 bash run.sh forward"
        exit 1
    fi
fi

# ── validate inputs ───────────────────────────────────────────────────────────
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory '$INPUT_DIR' not found."
    exit 1
fi

if ! compgen -G "$INPUT_DIR/inparam*.yaml" > /dev/null; then
    echo "ERROR: no inparam YAML files were found in $INPUT_DIR."
    exit 1
fi

# ── locate axisem3d binary ────────────────────────────────────────────────────
if [[ -f "./axisem3d" ]]; then
    AXISEM3D_BIN="$(pwd)/axisem3d"
elif command -v axisem3d &>/dev/null; then
    AXISEM3D_BIN="$(command -v axisem3d)"
else
    echo "ERROR: axisem3d binary not found in this directory or on PATH."
    echo "       Copy the compiled binary here or put it on PATH."
    exit 1
fi

# ── auto-detect MPI launcher from binary's linked libmpi ─────────────────────
# Mirrors the logic in axikernels.core.kernels.objective_function so the same
# MPI runtime is used for both forward and backward runs.
MPIRUN="mpirun"
LD_PREPEND=""

if command -v ldd &>/dev/null; then
    while IFS= read -r line; do
        if [[ "$line" == *libmpi* ]] && [[ "$line" == *"=>"* ]]; then
            lib_path=$(echo "$line" | awk -F '=>' '{print $2}' | awk '{print $1}')
            if [[ -n "$lib_path" ]] && [[ "$lib_path" != "not" ]]; then
                conda_lib="$(dirname "$lib_path")"
                candidate="$(dirname "$conda_lib")/bin/mpirun"
                if [[ -f "$candidate" ]]; then
                    MPIRUN="$candidate"
                    LD_PREPEND="$conda_lib"
                fi
                break
            fi
        fi
    done < <(ldd "$AXISEM3D_BIN" 2>/dev/null)
fi

if [[ "$MPIRUN" == */* ]]; then
    if [[ ! -x "$MPIRUN" ]]; then
        echo "ERROR: detected MPI launcher is not executable: $MPIRUN"
        exit 1
    fi
elif ! command -v "$MPIRUN" &>/dev/null; then
    echo "ERROR: MPI launcher '$MPIRUN' was not found on PATH."
    echo "       If your system uses a different launcher, adjust PATH or edit run.sh."
    exit 1
fi

# ── set up simulation directory ───────────────────────────────────────────────
mkdir -p "${SIMU_DIR}/input"
cp -r "${INPUT_DIR}/"* "${SIMU_DIR}/input/"
cp "$AXISEM3D_BIN" "${SIMU_DIR}/axisem3d"

# ── run forward simulation ────────────────────────────────────────────────────
echo "Running forward simulation in ${SIMU_DIR}/ with ${NRANKS} MPI ranks..."
echo "  binary : $AXISEM3D_BIN"
echo "  mpirun : $MPIRUN"
[[ -n "$LD_PREPEND" ]] && echo "  LD_LIBRARY_PATH prepend: $LD_PREPEND"

cd "${SIMU_DIR}"
if [[ -n "$LD_PREPEND" ]]; then
    export LD_LIBRARY_PATH="${LD_PREPEND}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
"$MPIRUN" -np "$NRANKS" ./axisem3d
cd ..

echo ""
echo "Forward simulation complete."
echo "Element output is in ${SIMU_DIR}/output/elements/"
echo ""
echo "Next step:"
echo "  python compute_kernels.py"
echo ""
echo "(The backward simulation will be created and run automatically.)"

