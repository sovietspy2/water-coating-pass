#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

INPUT_PDB="$1"
INPUT_DIR="$(dirname "$INPUT_PDB")"
PDB_NAME="$(basename "${INPUT_PDB%.pdb}")"

"$SCRIPT_DIR"/tinker-pdb-converter.sh "$INPUT_PDB"

cp -f -- $SCRIPT_DIR/amber99.prm $INPUT_DIR/amber99.prm

echo "${INPUT_DIR}/${PDB_NAME}.xyz"

SECONDS=0
TIMEFILE="${PDB_NAME}-mm-process-time.txt"

write_runtime() {
  printf 'runtime for mm is %s seconds\n' "$SECONDS" > "$TIMEFILE"
}

trap write_runtime EXIT

# --- 1. SETUP ---
pdbxyz "$INPUT_PDB" <<EOF
ALL
amber99.prm
EOF

# --- 2. MINIMIZATION ---
minimize "${INPUT_DIR}/${PDB_NAME}.xyz" <<EOF
amber99.prm
0.01
EOF


# --- 3. BUILD RESTRICTION KEY ---
awk '$3 == "CA" {print "RESTRICT  " $2 "  10"}' "${INPUT_PDB}" > "${INPUT_DIR}/key.key"

cat >> "${INPUT_DIR}/key.key" <<EOF
PARAMETERS amber99.prm

RATTLE

vdw-cutoff 9.0
chg-cutoff 9.0
EOF


# --- 4. DYNAMICS ---

dynamic "${INPUT_DIR}/${PDB_NAME}".xyz_2 -k "$INPUT_DIR/key.key" <<EOF
100000
1.0
10
2
300
EOF


# --- 5. CONVERT BACK TO PDB ---
xyzpdb "${INPUT_DIR}/${PDB_NAME}.arc" <<EOF
amber99.prm
PDB
EOF

OUT="${INPUT_DIR}/filtered_renum.pdb" # Double check
cp -f -- "${INPUT_DIR}/${PDB_NAME}.pdb_2" "$OUT"

echo "tinker.sh done"