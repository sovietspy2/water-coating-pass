#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

ns_to_steps() {
  local NS="${1:?usage: ns_to_steps <nanoseconds> [dt_fs]}"
  local DT_FS="${2:-1.0}"

  awk -v NS="$NS" -v DT="$DT_FS" 'BEGIN {
    printf "%.0f\n", (NS * 1000000) / DT
  }'
}

ns_to_ps() {
  local NS="${1:?usage: ns_to_ps <nanoseconds>}"
  awk -v NS="$NS" 'BEGIN {
    printf "%.6f\n", NS * 1000.0
  }'
}

arc_lines_per_frame() {
  local ARC_FILE="${1:?usage: arc_lines_per_frame <file.arc>}"

  local N_ATOMS SECOND_FIRST
  N_ATOMS=$(awk 'NR==1 {print $1; exit}' "$ARC_FILE")
  SECOND_FIRST=$(awk 'NR==2 {print $1; exit}' "$ARC_FILE")

  case "$SECOND_FIRST" in
    *[!0-9]* ) echo $((N_ATOMS + 2)) ;;
    * )        echo $((N_ATOMS + 1)) ;;
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

# Tinker9-GPU support: export TINKER_GPU=1 to use GPU-accelerated commands.
# Only dynamic and minimize have GPU counterparts; file-format utilities
# (pdbxyz, arcedit, xyzpdb) are unchanged — Tinker9 uses identical file formats.
if [[ -n "${TINKER_GPU:-}" ]]; then
    MINIMIZE_CMD="minimize9"
    DYNAMIC_CMD="dynamic9"
    command -v minimize9 >/dev/null 2>&1 || { log "ERROR: minimize9 not found in PATH (TINKER_GPU is set)"; exit 1; }
    command -v dynamic9  >/dev/null 2>&1 || { log "ERROR: dynamic9 not found in PATH (TINKER_GPU is set)";  exit 1; }
    log "TINKER_GPU is set: using Tinker9-GPU commands (minimize9, dynamic9)"
else
    MINIMIZE_CMD="minimize"
    DYNAMIC_CMD="dynamic"
    log "TINKER_GPU is not set: using standard Tinker CPU commands (minimize, dynamic)"
fi

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
run_step "$MINIMIZE_CMD" "${INPUT_DIR}/${PDB_NAME}.xyz" <<EOF
amber99.prm
0.01
EOF

# 5) Build key file
# We must not add constraint for WATER O or H
log "Step 5: Building restriction key for protein heavy atoms only (excluding all water atoms)"

awk '
function trim(S) {
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", S)
  return S
}

function is_hydrogen(ATOM, ELEM) {
  ATOM = toupper(ATOM)
  ELEM = toupper(ELEM)

  # Prefer element column when present.
  if (ELEM != "") {
    return (ELEM == "H")
  }

  # Fallback to atom name if needed.
  return (ATOM ~ /^[0-9]*H/)
}

function is_water_residue(RESN) {
  RESN = toupper(RESN)
  return (RESN ~ /^(HOH|WAT|SOL|H2O|DOD|D2O|TIP|OH2|OD2)$/)
}

# Only protein/polymer ATOM records; ignore HETATM entirely.
# Also exclude any water-like residue names just to be safe.
 /^ATOM  / {
  SERIAL = trim(substr($0, 7, 5))
  ATOM   = trim(substr($0, 13, 4))
  RESN   = trim(substr($0, 18, 3))
  ELEM   = trim(substr($0, 77, 2))

  if (!is_water_residue(RESN) && !is_hydrogen(ATOM, ELEM)) {
    print "RESTRICT  " SERIAL "  200"
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

# Tinker9-GPU: two extra key.key settings required for GPU execution.
# 1) Beeman integrator (CPU default) is not implemented in Tinker9-GPU — use velocity Verlet.
# 2) REMOVE-INERTIA 0: disables angular-momentum removal (mdrestRemoveAngularMomentum_cu is
#    unimplemented for non-PBC systems in Tinker9-GPU — harmless to skip for short runs).
if [[ -n "${TINKER_GPU:-}" ]]; then
  printf 'INTEGRATOR VERLET\nREMOVE-INERTIA 0\n' >> "${INPUT_DIR}/key.key"
  log "TINKER_GPU: added INTEGRATOR VERLET + REMOVE-INERTIA 0 to key.key"
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
  SAVE_INTERVAL_PS="$(awk -v STEPS="$SAVE_EVERY_STEPS" -v DT="$DT_FS" 'BEGIN {
    printf "%.6f\n", (STEPS * DT) / 1000.0
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

run_step "$DYNAMIC_CMD" "${INPUT_DIR}/${PDB_NAME}.xyz_2" -k "${INPUT_DIR}/key.key" <<EOF
${N_STEPS}
${DT_FS}
${SAVE_INTERVAL_PS}
2
300
EOF

# 7) Choose final structure source
ARC_FILE="${INPUT_DIR}/${PDB_NAME}.arc"
FINAL_XYZ=""

if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
  if [[ ! -f "$ARC_FILE" ]]; then
    log "ERROR: ARC file not found: $ARC_FILE"
    exit 1
  fi

  TOTAL_LINES=$(wc -l < "$ARC_FILE" | tr -d ' ')
  LINES_PER_FRAME="$(arc_lines_per_frame "$ARC_FILE")"
  LAST_FRAME_NUMBER=$((TOTAL_LINES / LINES_PER_FRAME))
  log "Last frame index in .arc file is ${LAST_FRAME_NUMBER}"

  log "Step 7: Extracting last frame with arcedit"
  run_step arcedit "$ARC_FILE" <<EOF
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
OUTPUT_PDB="${INPUT_DIR}/next_step.pdb"
cp -f -- "$LAST_PDB" "$OUTPUT_PDB"
log "Copied final structure to next_step.pdb"

if [[ "$MOBYWAT_OUTPUT_ENABLED" != "true" ]]; then
  log "MobyWat output disabled; skipping trajectory PDB generation"
  log "tinker.sh completed successfully"
  exit 0
fi

# 10) Build a multi-model trajectory PDB from the full ARC history
log "Step 10: Building multi-model trajectory PDB from ARC"
run_step xyzpdb "$ARC_FILE" <<EOF
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

  log "Making sure Reference PDB is compatible with mobywat."
  run_step "${SCRIPT_DIR}/reference-pdb-preprocessor.py" ${SYSTEM_REF_PDB}

  log "Adding mobywat params to ${SYSTEM_REF_PDB}"
  run_step "${SCRIPT_DIR}/add-mobywat-analysis-params.sh" ${SYSTEM_REF_PDB}

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