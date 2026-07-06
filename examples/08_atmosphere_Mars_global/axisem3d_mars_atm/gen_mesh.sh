#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PERIOD="${1:-20}"
OUTPUT_FILE="${2:-${SCRIPT_DIR}/input/tayak_60km_${PERIOD}s.e}"
PREP_BM="${SCRIPT_DIR}/tayak_60km_atm_meshing.bm"

python3 "${SCRIPT_DIR}/prepare_positive_atmosphere_bm.py" \
  --input "${SCRIPT_DIR}/tayak_60km_atm.bm" \
  --output "${PREP_BM}" \
  --max-atmosphere-step 100.0

python -m salvus_mesh_lite.interface AxiSEM \
  --basic.model "${PREP_BM}" \
  --basic.period "${PERIOD}" \
  --output_filename "${OUTPUT_FILE}"

python3 "${SCRIPT_DIR}/validate_exodus_fluid_rho.py" "${OUTPUT_FILE}"
