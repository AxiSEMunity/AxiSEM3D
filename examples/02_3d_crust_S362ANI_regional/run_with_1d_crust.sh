#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
	echo "Run this script instead of sourcing it." >&2
	return 1
fi

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
AXISEM3D_BIN="${AXISEM3D_BIN:-$REPO_DIR/build/axisem3d}"
RUN_DIR="$SCRIPT_DIR/simu_with_1d_crust"
INPUT_DIR="$RUN_DIR/input"
OUTPUT_DIR="${OUTPUT_DIR:-$RUN_DIR/output}"
NP="${NP:-4}"

if [[ -f /etc/profile.d/lmod.sh ]]; then
	set +u
	source /etc/profile.d/lmod.sh
	set -u
fi

if type module >/dev/null 2>&1; then
	module load ohpc
	module load gnu13/13.2.0
	module load openmpi5/5.0.5
	module load netcdf/4.9.2
	module load fftw/3.3.10
	module load metis/5.1.0
fi

if [[ ! -x "$AXISEM3D_BIN" ]]; then
	echo "AxiSEM3D binary not found at $AXISEM3D_BIN" >&2
	exit 1
fi

if ! command -v mpirun >/dev/null 2>&1; then
	echo "mpirun is not available; load the MPI runtime first." >&2
	exit 1
fi

mkdir -p "$RUN_DIR"
rm -rf "$INPUT_DIR"
mkdir -p "$INPUT_DIR"
cp -r "$SCRIPT_DIR/input_share/." "$INPUT_DIR/"
cp -r "$SCRIPT_DIR/input_with_1d_crust/." "$INPUT_DIR/"

mpirun -np "$NP" "$AXISEM3D_BIN" \
	--input "$INPUT_DIR" \
	--output "$OUTPUT_DIR"
