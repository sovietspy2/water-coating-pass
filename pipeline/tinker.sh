#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

INPUT_PDB="$(normalize_input_pdb "$1")"
INPUT_DIR="$(dirname "$INPUT_PDB")"
BASE="$(basename "${INPUT_PDB%.pdb}")"

SECONDS=0
TIMEFILE="${INPUT_PDB}-mm-process-time.txt"

write_runtime() {
  local elapsed="$SECONDS"
  # Write runtime (in seconds) to the requested file
  printf 'runtime for mm is %s seconds\n' "$elapsed" > "$TIMEFILE"
}
trap write_runtime EXIT

echo "tinker.sh start"

# --- 1. SETUP ---
pdbxyz "$INPUT_PDB" <<EOF
ALL
${SCRIPT_DIR}/amber99.prm
EOF

echo "pdbxyz done"

#echo "hello 2"

# --- 2. MINIMIZATION ---
minimize "${BASE}.xyz" <<EOF
${SCRIPT_DIR}/amber99.prm
0.01
EOF


# --- 3. BUILD RESTRICTION KEY ---
awk '$3 == "CA" {print "RESTRICT  " $2 "  10"}' "${INPUT_PDB}" > key.key

cat >> key.key <<EOF
PARAMETERS ${SCRIPT_DIR}/amber99.prm

RATTLE

vdw-cutoff 9.0
chg-cutoff 9.0
EOF

echo "awk cat done"

# --- 4. DYNAMICS ---
dynamic "${BASE}.xyz_2" -k key.key <<EOF
100000
1.0
10
2
300
EOF

echo "dynamic done"

# --- 5. CONVERT BACK TO PDB ---
xyzpdb "${BASE}.arc" <<EOF
${SCRIPT_DIR}/amber99.prm
PDB
EOF

echo "xyzpdb done"

OUT="${BASE}_mm.pdb"
cp -f -- "${BASE}.pdb_2" "$OUT"

echo "tinker.sh done"