#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  echo "Run this script instead of sourcing it." >&2
  return 1
fi

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
AXISEM3D_BIN="${AXISEM3D_BIN:-$REPO_DIR/build/axisem3d}"
CASE_REL="${1:-}"

if [[ -z "$CASE_REL" ]]; then
  echo "Usage: $0 <case-relative-path>" >&2
  exit 2
fi

if [[ ! -d "$SCRIPT_DIR/$CASE_REL" ]]; then
  echo "Case directory not found: $CASE_REL" >&2
  exit 1
fi

CASE_DIR="$(cd -- "$SCRIPT_DIR/$CASE_REL" && pwd)"
INPUT_DIR="$CASE_DIR/input"
OUTPUT_DIR="${OUTPUT_DIR:-$CASE_DIR/output}"
NP="${NP:-4}"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory not found at $INPUT_DIR" >&2
  exit 1
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

if ! command -v mpirun >/dev/null 2>&1; then
  echo "mpirun is not available; load the MPI runtime first." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

mpirun -np "$NP" "$AXISEM3D_BIN" \
  --input "$INPUT_DIR" \
  --output "$OUTPUT_DIR"