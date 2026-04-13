#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  echo "Run this script instead of sourcing it." >&2
  return 1
fi

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
AXISEM3D_BIN="${AXISEM3D_BIN:-$REPO_DIR/build/axisem3d}"
OUTPUT_DIR="${OUTPUT_DIR:-$SCRIPT_DIR/output}"
NP="${NP:-4}"

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

mpirun -np "$NP" "$AXISEM3D_BIN" \
  --input "$SCRIPT_DIR/input" \
  --output "$OUTPUT_DIR"