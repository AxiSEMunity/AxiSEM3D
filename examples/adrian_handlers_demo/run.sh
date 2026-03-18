#!/usr/bin/env bash
# Adrian Handlers Demo – run script
#
# Runs a tiny AxiSEM3D simulation and writes station + element output ready
# for post-processing with the axikernels handlers.
#
# Usage:
#   bash run.sh [--help] [--dry-run]
#
# Environment variables:
#   NRANKS   Number of MPI ranks (default: 2)
#
# Prerequisites:
#   - A compiled axisem3d binary placed in this directory or on PATH.
#   - An MPI implementation (mpirun / mpiexec) available in the environment.
#   - axikernels installed (needed for the notebook post-processing step only).
#
# The simulation writes output to sim_run/output/ (gitignored).
# Use clean.sh to remove that directory when done.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_DIR="${SCRIPT_DIR}/input"
RUN_DIR="${SCRIPT_DIR}/sim_run"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --help|-h)
            cat <<'EOF'
Usage:
  bash run.sh [--help] [--dry-run]

Options:
  --help      Show this message and exit.
  --dry-run   Validate input files and report what would be done without
              running the simulation.  Does not require the axisem3d binary.

Environment variables:
  NRANKS      Number of MPI ranks (default: 2)

Prerequisites:
  - axisem3d binary in this directory or on PATH.
  - MPI launcher (mpirun / mpiexec) available.
  - axikernels Python package installed (for the notebook step only).

Expected outputs (after a full run):
  sim_run/output/stations/GSN_Station_Grid/
  sim_run/output/elements/mantle/

Example:
  cp /path/to/build/axisem3d .
  NRANKS=4 bash run.sh
EOF
            exit 0
            ;;
        --dry-run)
            DRY_RUN=1
            ;;
        *)
            echo "Unknown argument: $arg  (use --help for usage)" >&2
            exit 1
            ;;
    esac
done

NRANKS="${NRANKS:-2}"
if ! [[ "${NRANKS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: NRANKS must be a positive integer (got: ${NRANKS})" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Validate input files
# ---------------------------------------------------------------------------
for required in \
    "${INPUT_DIR}/inparam.model.yaml" \
    "${INPUT_DIR}/inparam.output.yaml" \
    "${INPUT_DIR}/inparam.source.yaml" \
    "${INPUT_DIR}/inparam.advanced.yaml" \
    "${INPUT_DIR}/inparam.nr.yaml" \
    "${INPUT_DIR}/global_mesh__prem_ani__50s.e" \
    "${INPUT_DIR}/1dmodel_axisem.bm" \
    "${INPUT_DIR}/GSN_small.txt" \
    "${INPUT_DIR}/HANDLERS_EXAMPLE_cat.xml"; do
    if [[ ! -f "${required}" ]]; then
        echo "ERROR: required input file not found: ${required}" >&2
        exit 1
    fi
done
echo "Input files OK: ${INPUT_DIR}"

# ---------------------------------------------------------------------------
# Dry-run: report and exit without touching sim_run/ or checking the binary
# ---------------------------------------------------------------------------
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo ""
    echo "--- DRY RUN: all input files present, no simulation will be run ---"
    echo ""
    echo "Run directory  : ${RUN_DIR}"
    echo "MPI ranks      : ${NRANKS}"
    echo ""
    echo "NOTE: the axisem3d binary is NOT validated in dry-run mode."
    echo "      Copy it here or put it on PATH before the full run."
    echo ""
    echo "Expected outputs after a full run:"
    echo "  ${RUN_DIR}/output/stations/GSN_Station_Grid/"
    echo "  ${RUN_DIR}/output/elements/mantle/"
    exit 0
fi

# ---------------------------------------------------------------------------
# Locate axisem3d binary (user-provided; NOT bundled in this example)
# ---------------------------------------------------------------------------
if [[ -f "${SCRIPT_DIR}/axisem3d" ]]; then
    AXISEM3D_BIN="${SCRIPT_DIR}/axisem3d"
elif command -v axisem3d &>/dev/null; then
    AXISEM3D_BIN="$(command -v axisem3d)"
else
    echo "" >&2
    echo "ERROR: axisem3d binary not found." >&2
    echo "" >&2
    echo "  This example does NOT ship a pre-built binary." >&2
    echo "  Build AxiSEM3D following the main README, then either:" >&2
    echo "    cp /path/to/build/axisem3d ." >&2
    echo "  or add it to PATH." >&2
    exit 1
fi
echo "Binary : ${AXISEM3D_BIN}"

# ---------------------------------------------------------------------------
# Resolve MPI launcher
# ---------------------------------------------------------------------------
MPI_RUNNER=""
if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/mpirun" ]]; then
    MPI_RUNNER="${CONDA_PREFIX}/bin/mpirun"
elif [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/mpiexec" ]]; then
    MPI_RUNNER="${CONDA_PREFIX}/bin/mpiexec"
elif command -v mpirun &>/dev/null; then
    MPI_RUNNER="$(command -v mpirun)"
elif command -v mpiexec &>/dev/null; then
    MPI_RUNNER="$(command -v mpiexec)"
else
    echo "ERROR: no MPI launcher found. Activate the matching conda environment or install MPI." >&2
    exit 1
fi
echo "MPI    : ${MPI_RUNNER}"

# ---------------------------------------------------------------------------
# Set up a fresh run directory
# ---------------------------------------------------------------------------
echo ""
echo "Setting up run directory: ${RUN_DIR}"
mkdir -p "${RUN_DIR}"

# Remove only known AxiSEM3D output/log directories to avoid wiping anything
# the user may have placed alongside them.
for stale in output/ logs/; do
    if [[ -d "${RUN_DIR}/${stale}" ]]; then
        echo "  Removing stale: ${RUN_DIR}/${stale}"
        rm -rf "${RUN_DIR:?}/${stale}"
    fi
done

# Deploy binary and input into the run directory
cp "${AXISEM3D_BIN}" "${RUN_DIR}/axisem3d"
chmod +x "${RUN_DIR}/axisem3d"

rm -rf "${RUN_DIR}/input"
cp -r "${INPUT_DIR}" "${RUN_DIR}/input"

echo "Deployed to run directory."

# ---------------------------------------------------------------------------
# Run the simulation
# ---------------------------------------------------------------------------
echo ""
echo "Running AxiSEM3D with ${NRANKS} MPI rank(s) ..."
echo "  Working directory: ${RUN_DIR}"
echo "  MPI launcher     : ${MPI_RUNNER}"
echo ""

pushd "${RUN_DIR}" > /dev/null
"${MPI_RUNNER}" -n "${NRANKS}" ./axisem3d input
popd > /dev/null

echo ""
echo "Simulation complete."
echo "Outputs:"
echo "  ${RUN_DIR}/output/stations/GSN_Station_Grid/"
echo "  ${RUN_DIR}/output/elements/mantle/"
echo ""
echo "Next step: open handlers_demo.ipynb in Jupyter."
