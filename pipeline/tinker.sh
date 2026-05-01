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

# Setup logging in the INPUT_DIR (output directory)
LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging $LOGFILE

log "Starting tinker.sh (TINKER mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_DIR=$INPUT_DIR"
log "PDB_NAME=$PDB_NAME"

SECONDS=0
TIMEFILE="${INPUT_DIR}/${PDB_NAME}-mm-process-time.txt"

write_runtime() {
  printf 'runtime for mm is %s seconds\n' "$SECONDS" > "$TIMEFILE"
  log "Runtime written to $TIMEFILE: $SECONDS seconds"
}

trap write_runtime EXIT

touch "${INPUT_DIR}/TINKER.protocol"
log "Created TINKER.protocol file"

# Reformatting PDB so tinker can use it
log "Step 1: Changing pdb format to fit tinker requirements"
run_step "$SCRIPT_DIR"/format-pdb.sh "$INPUT_PDB"

log "Step 2: Copying amber99.prm to input directory"
cp -f -- $SCRIPT_DIR/amber99.prm $INPUT_DIR/amber99.prm
log "Copied amber99.prm to $INPUT_DIR/"

log "Generated XYZ file: ${INPUT_DIR}/${PDB_NAME}.xyz"

# --- 1. SETUP --- A
log "Step 3: Running pdbxyz (PDB to XYZ conversion) A VERSION"
run_step pdbxyz "$INPUT_PDB" <<EOF
ALL
amber99.prm
EOF

# --- 2. MINIMIZATION ---
log "Step 4: Running minimize (energy minimization)"
run_step minimize "${INPUT_DIR}/${PDB_NAME}.xyz" <<EOF
amber99.prm
0.01
EOF

log "Step 5: Building restriction key for CA atoms"

# --- 3. BUILD RESTRICTION KEY ---
awk '$3 == "CA" {print "RESTRICT  " $2 "  10"}' "${INPUT_PDB}" > "${INPUT_DIR}/key.key"
log "Generated key.key with CA restrictions"

cat >> "${INPUT_DIR}/key.key" <<EOF
PARAMETERS amber99.prm

RATTLE

vdw-cutoff 9.0
chg-cutoff 9.0
EOF

log "Appended PARAMETERS, RATTLE, and cutoff settings to key.key"

# --- 4. DYNAMICS ---
log "Step 6: Running dynamic (molecular dynamics simulation)"
run_step dynamic "${INPUT_DIR}/${PDB_NAME}".xyz_2 -k "$INPUT_DIR/key.key" <<EOF
100000
1.0
10
2
300
EOF

# Extracting last frame of .arc file
NATOMS=$(awk 'NR==1 {print $1; exit}' "$ARC")
SECOND_FIRST=$(awk 'NR==2 {print $1; exit}' "$ARC")
TOTAL_LINES=$(wc -l < "$ARC" | tr -d ' ')

case "$SECOND_FIRST" in
    *[!0-9]* ) LINES_PER_FRAME=$((NATOMS + 2)) ;;
    * )        LINES_PER_FRAME=$((NATOMS + 1)) ;;
esac

LAST_FRAME_NUMBER=$((TOTAL_LINES / LINES_PER_FRAME))
log "Last frame in .arc file is ${LAST_FRAME}"

run_step arcedit"${INPUT_DIR}/${PDB_NAME}.arc" <<EOF
2
${LAST_FRAME_NUMBER} ${LAST_FRAME_NUMBER} 1

EOF

LAST_FRAME="${INPUT_DIR}/${PDB_NAME}.0${LAST_FRAME_NUMBER}"
log "Expected single frame arc file is: ${}"

# --- 5. CONVERT BACK TO PDB ---
log "Step 7: Running xyzpdb (converting trajectory arc to PDB)"
run_step xyzpdb ${LAST_FRAME} <<EOF
amber99.prm
PDB
EOF

log "Step 8: Preparing final output"
OUT="${INPUT_DIR}/filtered_renum.pdb"
cp -f -- "${INPUT_DIR}/${PDB_NAME}.pdb_2" "$OUT"
log "Copied final structure to filtered_renum.pdb"

log "Step 9: Post processing PDB"
run_step "$SCRIPT_DIR"/format-pdb.sh "${INPUT_DIR}/filtered_renum.pdb"

log "tinker.sh completed successfully"