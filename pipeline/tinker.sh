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

ns_to_ps() {
  local ns="${1:?usage: ns_to_ps <nanoseconds>}"
  awk -v ns="$ns" 'BEGIN {
    printf "%.6f\n", ns * 1000.0
  }'
}

arc_lines_per_frame() {
  local arc="${1:?usage: arc_lines_per_frame <file.arc>}"

  local natoms second_first
  natoms=$(awk 'NR==1 {print $1; exit}' "$arc")
  second_first=$(awk 'NR==2 {print $1; exit}' "$arc")

  case "$second_first" in
    *[!0-9]* ) echo $((natoms + 2)) ;;
    * )        echo $((natoms + 1)) ;;
  esac
}

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <input.pdb> <md_duration_ns> <mobywat_output_enabled:true|false> <target_frames> <reference.pdb>" >&2
  exit 1
fi

INPUT_PDB="$1"
INPUT_DIR="$(dirname "$INPUT_PDB")"
PDB_NAME="$(basename "${INPUT_PDB%.pdb}")"
MD_DURATION="$2"
MOBYWAT_OUTPUT_ENABLED="$3"
TARGET_FRAMES="${4:-1000}"
REFERENCE_PDB="${5:-}"
DT_FS="1.0" # 1 femtosecond

LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging "$LOGFILE"

log "Starting tinker.sh (TINKER mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_DIR=$INPUT_DIR"
log "PDB_NAME=$PDB_NAME"
log "MD_DURATION=$MD_DURATION"
log "MOBYWAT_OUTPUT_ENABLED=$MOBYWAT_OUTPUT_ENABLED"
log "TARGET_FRAMES=$TARGET_FRAMES"
log "DT_FS=$DT_FS"
log "REFERENCE_PDB=$REFERENCE_PDB"

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
run_step "$SCRIPT_DIR/format-pdb.sh" "$INPUT_PDB"

# 2) Parameter file
log "Step 2: Copying amber99.prm to input directory"
cp -f -- "$SCRIPT_DIR/amber99.prm" "$INPUT_DIR/amber99.prm"
log "Copied amber99.prm to $INPUT_DIR/"

# 3) PDB -> XYZ
log "Step 3: Running pdbxyz (PDB to XYZ conversion)"
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

# 5) Build key file
# We must not add constraint for WATER O or H
log "Step 5: Building restriction key for protein heavy atoms only (excluding all water atoms)"

awk '
function trim(s) {
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
  return s
}

function is_hydrogen(atom, elem) {
  atom = toupper(atom)
  elem = toupper(elem)

  # Prefer element column when present.
  if (elem != "") {
    return (elem == "H")
  }

  # Fallback to atom name if needed.
  return (atom ~ /^[0-9]*H/)
}

function is_water_residue(resn) {
  resn = toupper(resn)
  return (resn ~ /^(HOH|WAT|SOL|H2O|DOD|D2O|TIP|OH2|OD2)$/)
}

# Only protein/polymer ATOM records; ignore HETATM entirely.
# Also exclude any water-like residue names just to be safe.
 /^ATOM  / {
  serial = trim(substr($0, 7, 5))
  atom   = trim(substr($0, 13, 4))
  resn   = trim(substr($0, 18, 3))
  elem   = trim(substr($0, 77, 2))

  if (!is_water_residue(resn) && !is_hydrogen(atom, elem)) {
    print "RESTRICT  " serial "  200"
  }
}
' "$INPUT_PDB" > "${INPUT_DIR}/key.key"

log "Generated key.key with restrictions for protein heavy atoms only"

cat >> "${INPUT_DIR}/key.key" <<EOF
PARAMETERS amber99.prm

RATTLE
vdw-cutoff 9.0
chg-cutoff 9.0
EOF

if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
  echo "ARCHIVE" >> "${INPUT_DIR}/key.key"
  log "Trajectory mode enabled: ARCHIVE written to key.key"
else
  echo "NO-ARCHIVE" >> "${INPUT_DIR}/key.key"
  log "Fast mode enabled: NO-ARCHIVE written to key.key"
fi

# 6) Dynamics setup
N_STEPS="$(ns_to_steps "$MD_DURATION" "$DT_FS")"
TOTAL_TIME_PS="$(ns_to_ps "$MD_DURATION")"

if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
  SAVE_EVERY_STEPS=$(( N_STEPS / TARGET_FRAMES ))
  if (( SAVE_EVERY_STEPS < 1 )); then
    SAVE_EVERY_STEPS=1
  fi

  ACTUAL_FRAMES=$(( (N_STEPS + SAVE_EVERY_STEPS - 1) / SAVE_EVERY_STEPS ))
  SAVE_INTERVAL_PS="$(awk -v steps="$SAVE_EVERY_STEPS" -v dt="$DT_FS" 'BEGIN {
    printf "%.6f\n", (steps * dt) / 1000.0
  }')"
else
  SAVE_EVERY_STEPS="$N_STEPS"
  ACTUAL_FRAMES=1
  SAVE_INTERVAL_PS="$TOTAL_TIME_PS"
fi

log "Step 6: Running dynamic (molecular dynamics simulation)"
log "MD_DURATION=$MD_DURATION ns"
log "TOTAL_TIME_PS=$TOTAL_TIME_PS ps"
log "N_STEPS=$N_STEPS"
log "SAVE_EVERY_STEPS=$SAVE_EVERY_STEPS steps"
log "SAVE_INTERVAL_PS=$SAVE_INTERVAL_PS ps"
log "Expected saved frames ~= $ACTUAL_FRAMES"

run_step dynamic "${INPUT_DIR}/${PDB_NAME}.xyz_2" -k "${INPUT_DIR}/key.key" <<EOF
${N_STEPS}
${DT_FS}
${SAVE_INTERVAL_PS}
2
300
EOF

# 7) Choose final structure source
ARC="${INPUT_DIR}/${PDB_NAME}.arc"
FINAL_XYZ=""

if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
  if [[ ! -f "$ARC" ]]; then
    log "ERROR: ARC file not found: $ARC"
    exit 1
  fi

  TOTAL_LINES=$(wc -l < "$ARC" | tr -d ' ')
  LINES_PER_FRAME="$(arc_lines_per_frame "$ARC")"
  LAST_FRAME_NUMBER=$((TOTAL_LINES / LINES_PER_FRAME))
  log "Last frame index in .arc file is ${LAST_FRAME_NUMBER}"

  log "Step 7: Extracting last frame with arcedit"
  run_step arcedit "$ARC" <<EOF
2
${LAST_FRAME_NUMBER} ${LAST_FRAME_NUMBER} 1

EOF

  FINAL_XYZ="${INPUT_DIR}/${PDB_NAME}.$(printf '%03d' "$LAST_FRAME_NUMBER")"
  log "Expected single-frame XYZ is: ${FINAL_XYZ}"

  if [[ ! -f "$FINAL_XYZ" ]]; then
    log "WARNING: Expected XYZ file ${FINAL_XYZ} not found; listing directory:"
    ls -1 "$INPUT_DIR"
    log "Please adjust the naming pattern for arcedit output."
    exit 1
  fi
else
  # In fast mode, dynamic updates the input coordinates; use the latest XYZ directly.
  if [[ -f "${INPUT_DIR}/${PDB_NAME}.xyz_2" ]]; then
    FINAL_XYZ="${INPUT_DIR}/${PDB_NAME}.xyz_2"
  elif [[ -f "${INPUT_DIR}/${PDB_NAME}.xyz" ]]; then
    FINAL_XYZ="${INPUT_DIR}/${PDB_NAME}.xyz"
  else
    log "ERROR: No final XYZ file found after dynamics"
    exit 1
  fi
  log "Fast mode final XYZ selected: ${FINAL_XYZ}"
fi

# 8) XYZ -> PDB for final structure
log "Step 8: Running xyzpdb (converting final-frame XYZ to PDB)"
run_step xyzpdb "$FINAL_XYZ" <<EOF
amber99.prm
PDB
EOF

LAST_PDB=$(ls -1t "${INPUT_DIR}/${PDB_NAME}.pdb"* 2>/dev/null | head -n 1 || true)
if [[ -z "$LAST_PDB" ]]; then
  log "ERROR: No PDB file produced by xyzpdb"
  exit 1
fi
log "Selected final PDB: ${LAST_PDB}"

log "Step 9: Preparing final output"
OUT="${INPUT_DIR}/next_step.pdb"
cp -f -- "$LAST_PDB" "$OUT"
log "Copied final structure to next_step.pdb"

if [[ "$MOBYWAT_OUTPUT_ENABLED" != "true" ]]; then
  log "MobyWat output disabled; skipping trajectory PDB generation"
  log "tinker.sh completed successfully"
  exit 0
fi

# 10) Build a multi-model trajectory PDB from the full ARC history
log "Step 10: Building multi-model trajectory PDB from ARC"
run_step xyzpdb "$ARC" <<EOF
amber99.prm
PDB
EOF

log "Step 11: Reformatting PDB file"
# xyzpdb produces NAME.pdb_3
MOBYWAT_INPUT_PDB="system_mdl.pdb"
cp -f -- "${INPUT_DIR}/${PDB_NAME}.pdb_3" "${INPUT_DIR}/${MOBYWAT_INPUT_PDB}"


run_step "$SCRIPT_DIR/format-pdb.sh" "${INPUT_DIR}/${MOBYWAT_INPUT_PDB}"
log "${INPUT_DIR}/${MOBYWAT_INPUT_PDB} is ready to be processed by MobyWat!"

log "Step 12: Running mobywat!"
if [[ -n "${REFERENCE_PDB:-}" ]]; then
  log "REFERENCE_PDB is present and non-empty: $REFERENCE_PDB, VALIDATION MODE!"

  SYSTEM_REF_PDB="system_ref.pdb"
  log "${REFERENCE_PDB} is copied to: ${SYSTEM_REF_PDB}, will edit it!"
  cp ${REFERENCE_PDB} ${SYSTEM_REF_PDB} #creating local version, we will modify it

  log "Adding mobywat params to ${SYSTEM_REF_PDB}"
  run_step "${SCRIPT_DIR}"/add-mobywat-analysis-params.sh ${SYSTEM_REF_PDB}

  log "Making sure Reference PDB is compatible with mobywat."
  run_step "${SCRIPT_DIR}/reference-pdb-preprocessor.py" ${SYSTEM_REF_PDB}

  log "Remove TER operator ID if present from ${SYSTEM_REF_PDB}"
  run_step "${SCRIPT_DIR}"/remove-ter-id.sh ${SYSTEM_REF_PDB}

  log "Running mobywat validation"
  run_step mobywat -f ${MOBYWAT_INPUT_PDB} -r ${SYSTEM_REF_PDB} -t [A] -w Auto -n 1-1000 -cls MER -m Analysis -v Diagnostic
else
  log "REFERENCE_PDB is missing or empty, PREDICTION MODE!"
  log "Running mobywat prediction"
  run_step mobywat -f ${MOBYWAT_INPUT_PDB} -t [A] -w Auto -n 1-1000 -cls MER -m Prediction -v Diagnostic
fi

log "tinker.sh completed successfully"