#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  echo "Run this script instead of sourcing it." >&2
  return 1
fi

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(cd -- "$TOOLS_DIR/../.." && pwd)"
VENV_DIR="$TOOLS_DIR/.venv"
AXISEM3D_BIN="${AXISEM3D_BIN:-$REPO_DIR/build/axisem3d}"
OUTPUT_DIR="${OUTPUT_DIR:-$SCRIPT_DIR/output}"
NP="${NP:-16}"

if [[ ! -d "$VENV_DIR" ]]; then
  echo "Python environment not found at $VENV_DIR" >&2
  exit 1
fi

source "$VENV_DIR/bin/activate"

if [[ ! -f "$SCRIPT_DIR/input/paper_example.e" ]]; then
  (cd "$SCRIPT_DIR/input" && bash gen_mesh.sh)
fi

if [[ ! -f "$SCRIPT_DIR/input/paper_example.nc" ]]; then
  (cd "$SCRIPT_DIR" && python -u generate_3D_model.py)
fi

if [[ -f /etc/profile.d/lmod.sh ]]; then
  set +u
  source /etc/profile.d/lmod.sh
  set -u
fi

if type module >/dev/null 2>&1; then
  # Adapt these module names to your cluster environment.
  # The solver needs an MPI runtime, NetCDF-C, FFTW3, and METIS.
  # module load compiler mpi netcdf fftw metis
  :
fi

if [[ ! -x "$AXISEM3D_BIN" ]]; then
  echo "AxiSEM3D binary not found at $AXISEM3D_BIN" >&2
  exit 1
fi

mpirun -np "$NP" "$AXISEM3D_BIN" \
  --input "$SCRIPT_DIR/input" \
  --output "$OUTPUT_DIR"