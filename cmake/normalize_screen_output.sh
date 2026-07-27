#!/usr/bin/env bash

set -euo pipefail

input_file="$1"
output_file="$2"

grep -v "wall-clock time" "$input_file" | grep -v "wave propagation time" | grep -v " w-c time" | grep -v "rank-to-rank communications" > "$output_file"
