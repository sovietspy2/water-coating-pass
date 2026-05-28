#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

ns_to_steps() {
  local ns="${1:?usage: ns_to_steps <nanoseconds> [dt_fs]}"
  local dt_fs="${2:-1.0}"

  awk -v ns="$ns" -v dt="$dt_fs" 'BEGIN {
    printf "%.0f\n", (ns * 1000000) / dt
  }'
}

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

INPUT_PDB="$1"
INPUT_DIR="$(dirname "$INPUT_PDB")"
PDB_NAME="$(basename "${INPUT_PDB%.pdb}")"
MD_DURATION="$2"
MOBYWAT_OUTPUT_ENABLED="$3"
TARGET_FRAMES="${4:-100}" # This defines how many frames we want to save during simulation
DT_FS="1.0" # tinker dynamic related property

# Setup logging in the INPUT_DIR (output directory)
LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging "$LOGFILE"

log "Starting tinker.sh (TINKER mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_DIR=$INPUT_DIR"
log "PDB_NAME=$PDB_NAME"
log "MD_DURATION=$MD_DURATION"
log "MOBYWAT_OUTPUT_ENABLED=$MOBYWAT_OUTPUT_ENABLED"
log "TARGET_FRAMES=$TARGET_FRAMES"

SECONDS=0
TIMEFILE="${INPUT_DIR}/${PDB_NAME}-mm-process-time.txt"

write_runtime() {
  printf 'runtime for mm is %s seconds\n' "$SECONDS" > "$TIMEFILE"
  log "Runtime written to $TIMEFILE: $SECONDS seconds"
}

trap write_runtime EXIT

touch "${INPUT_DIR}/TINKER.protocol"
log "Created TINKER.protocol file"

# 1) PDB pre-formatting
log "Step 1: Changing pdb format to fit tinker requirements"
run_step "$SCRIPT_DIR"/format-pdb.sh "$INPUT_PDB"

# 2) Parameter file
log "Step 2: Copying amber99.prm to input directory"
cp -f -- $SCRIPT_DIR/amber99.prm $INPUT_DIR/amber99.prm
log "Copied amber99.prm to $INPUT_DIR/"

# --- 1. SETUP --- A
log "Step 3: Running pdbxyz (PDB to XYZ conversion) A VERSION"
run_step pdbxyz "$INPUT_PDB" <<EOF
ALL
amber99.prm
EOF
log "Generated XYZ file: ${INPUT_DIR}/${PDB_NAME}.xyz"

# 4) Minimization
log "Step 4: Running minimize (energy minimization)"
run_step minimize "${INPUT_DIR}/${PDB_NAME}.xyz" <<EOF
amber99.prm
0.01
EOF

# 5) CA-restrained key
log "Step 5: Building restriction key for CA atoms"
awk '$3 == "CA" {print "RESTRICT  " $2 "  10"}' "${INPUT_PDB}" > "${INPUT_DIR}/key.key"
log "Generated key.key with CA restrictions"

cat >> "${INPUT_DIR}/key.key" <<EOF
PARAMETERS amber99.prm

RATTLE
ARCHIVE

vdw-cutoff 9.0
chg-cutoff 9.0
EOF

log "Appended PARAMETERS, RATTLE, and cutoff settings to key.key"

# 6) Dynamics

# Convertion ns to steps
N_STEPS="$(ns_to_steps "$MD_DURATION" "$DT_FS")"

# Choose approximate number of saved frames in STEP units first
SAVE_EVERY_STEPS=$(( N_STEPS / TARGET_FRAMES ))
if (( SAVE_EVERY_STEPS < 1 )); then
  SAVE_EVERY_STEPS=1
fi

ACTUAL_FRAMES=$(( (N_STEPS + SAVE_EVERY_STEPS - 1) / SAVE_EVERY_STEPS ))

# TINKER dynamic expects TIME BETWEEN SAVES IN PICOSECONDS
SAVE_INTERVAL_PS="$(awk -v steps="$SAVE_EVERY_STEPS" -v dt="$DT_FS" 'BEGIN {
  printf "%.6f\n", (steps * dt) / 1000.0
}')"

log "Step 6: Running dynamic (molecular dynamics simulation)"
log "MD_DURATION=$MD_DURATION ns"
log "N_STEPS=$N_STEPS"
log "TARGET_FRAMES=$TARGET_FRAMES"
log "SAVE_EVERY_STEPS=$SAVE_EVERY_STEPS steps"
log "Expected saved frames ~= $ACTUAL_FRAMES"
log "SAVE_INTERNAL_PS=$SAVE_INTERVAL_PS"

run_step dynamic "${INPUT_DIR}/${PDB_NAME}.xyz_2" -k "$INPUT_DIR/key.key" <<EOF
${N_STEPS}
${DT_FS}
${SAVE_INTERVAL_PS}
2
300
EOF

# 7) Find last frame number in .arc
ARC="${INPUT_DIR}/${PDB_NAME}.arc"
if [[ ! -f "$ARC" ]]; then
  log "ERROR: ARC file not found: $ARC"
  exit 1
fi

NATOMS=$(awk 'NR==1 {print $1; exit}' "$ARC")
SECOND_FIRST=$(awk 'NR==2 {print $1; exit}' "$ARC")
TOTAL_LINES=$(wc -l < "$ARC" | tr -d ' ')

case "$SECOND_FIRST" in
    *[!0-9]* ) LINES_PER_FRAME=$((NATOMS + 2)) ;;
    * )        LINES_PER_FRAME=$((NATOMS + 1)) ;;
esac

LAST_FRAME_NUMBER=$((TOTAL_LINES / LINES_PER_FRAME))
log "Last frame index in .arc file is ${LAST_FRAME_NUMBER}"
## TODO: We need to refactor this, keeping all the frames for last run!!!

# 8) Extract last frame with ARCEDIT
log "Step 7: Extracting last frame with arcedit"
run_step arcedit "$ARC" <<EOF
2
${LAST_FRAME_NUMBER} ${LAST_FRAME_NUMBER} 1

EOF

# Guess output XYZ name; TINKER typically uses basename.NNN
LAST_FRAME_XYZ="${INPUT_DIR}/${PDB_NAME}.$(printf '%03d' "$LAST_FRAME_NUMBER")"
log "Expected single-frame XYZ is: ${LAST_FRAME_XYZ}"

if [[ ! -f "$LAST_FRAME_XYZ" ]]; then
  log "WARNING: Expected XYZ file ${LAST_FRAME_XYZ} not found; listing directory:"
  ls -1 "${INPUT_DIR}"
  log "Please adjust the naming pattern for arcedit output."
  exit 1
fi

# 9) XYZ -> PDB for that last frame
log "Step 8: Running xyzpdb (converting last-frame XYZ to PDB)"
run_step xyzpdb "$LAST_FRAME_XYZ" <<EOF
amber99.prm
PDB
EOF

# Find the newest PDB produced by xyzpdb
LAST_PDB=$(ls -1t "${INPUT_DIR}/${PDB_NAME}.pdb"* 2>/dev/null | head -n 1 || true)
if [[ -z "${LAST_PDB}" ]]; then
  log "ERROR: No PDB file produced by xyzpdb"
  exit 1
fi
log "Selected final PDB: ${LAST_PDB}"

log "Step 9: Preparing final output"
OUT="${INPUT_DIR}/filtered_renum.pdb"
cp -f -- "$LAST_PDB" "$OUT"
log "Copied final structure to filtered_renum.pdb"

if [[ "$MOBYWAT_OUTPUT_ENABLED" != "true" ]]; then
  log "Not building MobyWat output yet."
  exit 0
fi

# 10) Build a multi-model trajectory PDB from the full ARC history (OPTIONAL)
log "Step 10: Building multi-model trajectory PDB from ARC"
run_step "$SCRIPT_DIR/arc_to_pdb.sh" $ARC mobywat_input.pdb

log "tinker.sh completed successfully"
