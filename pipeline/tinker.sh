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

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.pdb> <md_duration_ns> <target_frames> <reference.pdb> <run_mobywat>" >&2
  echo "  md_duration_ns > 0 runs MD; 0 skips it (minimize only)." >&2
  echo "  run_mobywat (default 1): 1 runs MobyWat after MD; 0 skips it (intermediate cycle)." >&2
  exit 1
fi

INPUT_PDB="$1"
INPUT_DIR="$(dirname "$INPUT_PDB")"
PDB_NAME="$(basename "${INPUT_PDB%.pdb}")"
MD_DURATION="$2"
TARGET_FRAMES="${3:-1000}"
REFERENCE_PDB="${4:-}"
RUN_MOBYWAT="${5:-1}" # 1 = run MobyWat after MD (final cycle); 0 = skip (intermediate cycle)
DT_FS="1.0" # 1 femtosecond

LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging "$LOGFILE"

log "Starting tinker.sh (TINKER mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_DIR=$INPUT_DIR"
log "PDB_NAME=$PDB_NAME"
log "MD_DURATION=$MD_DURATION"
log "TARGET_FRAMES=$TARGET_FRAMES"
log "RUN_MOBYWAT=$RUN_MOBYWAT"
log "DT_FS=$DT_FS"
log "REFERENCE_PDB=$REFERENCE_PDB"

# Tinker9-GPU support: export TINKER_GPU=1 to use GPU-accelerated commands.
# Only dynamic and minimize have GPU counterparts; file-format utilities
# (pdbxyz, arcedit, xyzpdb) are unchanged — Tinker9 uses identical file formats.
case "${TINKER_GPU:-0}" in
    ''|0|false|no|off) TINKER_GPU_ENABLED=false ;;
    *)                 TINKER_GPU_ENABLED=true ;;
esac

if [[ "$TINKER_GPU_ENABLED" == true ]]; then
    DYNAMIC_CMD="dynamic9"
    command -v dynamic9  >/dev/null 2>&1 || { log "ERROR: dynamic9 not found in PATH (TINKER_GPU is set)";  exit 1; }
    log "TINKER_GPU is set: using Tinker9-GPU commands (dynamic9)"
else
    DYNAMIC_CMD="dynamic"
    log "TINKER_GPU is not set: using standard Tinker CPU commands (dynamic)"
fi

case "${MOBYWAT_DEBUG:-0}" in
    ''|0|false|no|off) MOBYWAT_DEBUG_ENABLED=false ;;
    *)                 MOBYWAT_DEBUG_ENABLED=true ;;
esac

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
run_step "$SCRIPT_DIR/format_pdb.py" "$INPUT_PDB"

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

log "Step 4: Building key file with protein restraints (before minimization)"

# Count protein atoms: all XYZ atoms before the first water (OW/HW).
PROTEIN_ATOMS=$(awk 'NR>1 && ($2=="OW" || $2=="HW") {print NR-2; exit}' \
  "${INPUT_DIR}/${PDB_NAME}.xyz")

if [[ -z "$PROTEIN_ATOMS" || "$PROTEIN_ATOMS" -lt 1 ]]; then
  log "ERROR: could not determine protein atom count from XYZ file"
  exit 1
fi

log "Protein atom count (heavy+H, before first water in XYZ): ${PROTEIN_ATOMS}"

cat > "${INPUT_DIR}/minimize.key" <<EOF
RESTRAIN-POSITION -1 ${PROTEIN_ATOMS} 2.0
PARAMETERS amber99.prm
OPENMP-THREADS 1
VDW-CUTOFF 9.0
CHARGE-CUTOFF 9.0
EOF

log "Generated minimize.key with: RESTRAIN-POSITION -1 ${PROTEIN_ATOMS} 2.0, OPENMP-THREADS 1, VDW-CUTOFF 9.0, CHARGE-CUTOFF 9.0"

# 5) Minimization — pass minimze.key so protein restraints are active.
#    minimize9 reads PARAMETERS from the key file; stdin only needs the gradient criterion.
log "Step 5: Running minimize (energy minimization with protein restraints)"
run_step minimize "${INPUT_DIR}/${PDB_NAME}.xyz" -k "${INPUT_DIR}/minimize.key" <<EOF
0.01
EOF

# 6) Dynamics + final structure — only when a positive MD duration is set.
# With duration 0 we skip dynamics entirely and take the minimized structure (.xyz_2) as the
# final frame (used by the intermediate cycles of a multi-iteration run: wdrop + minimize only).
ARC_FILE="${INPUT_DIR}/${PDB_NAME}.arc"
FINAL_XYZ=""

if md_enabled "$MD_DURATION"; then
  # Append dynamics-specific settings after minimization completes.
  cat >> "${INPUT_DIR}/md.key" <<EOF
PARAMETERS amber99.prm
RANDOMSEED 28480426
RESTRAIN-POSITION -1 ${PROTEIN_ATOMS} 300.0

RATTLE
VDW-CUTOFF 9.0
CHARGE-CUTOFF 9.0
EOF
  echo "ARCHIVE" >> "${INPUT_DIR}/md.key"
  log "Generated md.key with: RANDOMSEED 28480426, RESTRAIN-POSITION -1 ${PROTEIN_ATOMS} 300.0 (trajectory mode: ARCHIVE)"

  # Tinker9-GPU: two extra key settings required for GPU execution.
  # 1) Beeman integrator (CPU default) is not implemented in Tinker9-GPU — use velocity Verlet.
  # 2) REMOVE-INERTIA 0: disables angular-momentum removal (mdrestRemoveAngularMomentum_cu is
  #    unimplemented for non-PBC systems in Tinker9-GPU — harmless to skip for short runs).
  if [[ "$TINKER_GPU_ENABLED" == true ]]; then
    printf 'INTEGRATOR VERLET\nREMOVE-INERTIA 0\n' >> "${INPUT_DIR}/md.key"
    log "TINKER_GPU: added INTEGRATOR VERLET + REMOVE-INERTIA 0 to md.key"
  fi

  # Dynamics setup
  N_STEPS="$(ns_to_steps "$MD_DURATION" "$DT_FS")"
  TOTAL_TIME_PS="$(ns_to_ps "$MD_DURATION")"

  SAVE_EVERY_STEPS=$(( N_STEPS / TARGET_FRAMES ))
  if (( SAVE_EVERY_STEPS < 1 )); then
    SAVE_EVERY_STEPS=1
  fi

  ACTUAL_FRAMES=$(( (N_STEPS + SAVE_EVERY_STEPS - 1) / SAVE_EVERY_STEPS ))
  MOBYWAT_FRAMES=$(( N_STEPS / SAVE_EVERY_STEPS ))
  SAVE_INTERVAL_PS="$(awk -v STEPS="$SAVE_EVERY_STEPS" -v DT="$DT_FS" 'BEGIN {
    printf "%.6f\n", (STEPS * DT) / 1000.0
  }')"

  log "Step 6: Running dynamic (molecular dynamics simulation)"
  log "MD_DURATION=$MD_DURATION ns"
  log "TOTAL_TIME_PS=$TOTAL_TIME_PS ps"
  log "N_STEPS=$N_STEPS"
  log "SAVE_EVERY_STEPS=$SAVE_EVERY_STEPS steps"
  log "SAVE_INTERVAL_PS=$SAVE_INTERVAL_PS ps"
  log "Expected saved frames ~= $ACTUAL_FRAMES (MobyWat range 1-$MOBYWAT_FRAMES)"

  run_step "$DYNAMIC_CMD" "${INPUT_DIR}/${PDB_NAME}.xyz_2" -k "${INPUT_DIR}/md.key" <<EOF
${N_STEPS}
${DT_FS}
${SAVE_INTERVAL_PS}
2
300
EOF

  # 7) Choose final structure source: last frame of the trajectory
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
  # MD disabled: the minimized structure (.xyz_2, written by minimize) is the final frame.
  log "Step 6: MD disabled (MD_DURATION=0); using minimized structure as final frame"
  FINAL_XYZ="${INPUT_DIR}/${PDB_NAME}.xyz_2"
  if [[ ! -f "$FINAL_XYZ" ]]; then
    log "ERROR: minimized XYZ not found: $FINAL_XYZ"
    exit 1
  fi
  log "Final XYZ selected: ${FINAL_XYZ}"
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

if ! should_run_mobywat "$MD_DURATION" "$RUN_MOBYWAT"; then
  if ! md_enabled "$MD_DURATION"; then
    log "MD disabled (MD_DURATION=0); skipping trajectory PDB generation and MobyWat"
  else
    log "Intermediate cycle (RUN_MOBYWAT=0); MD ran but skipping trajectory PDB generation and MobyWat"
  fi
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

run_step "$SCRIPT_DIR/format_pdb.py" "${INPUT_DIR}/${MOBYWAT_INPUT_PDB}"
log "${INPUT_DIR}/${MOBYWAT_INPUT_PDB} is ready to be processed by MobyWat!"

log "Step 12: Running mobywat!"
MOBYWAT_ARGS=(-f "${MOBYWAT_INPUT_PDB}") #list of params to mobywat
if [[ -n "${REFERENCE_PDB:-}" ]]; then
  log "REFERENCE_PDB is present and non-empty: $REFERENCE_PDB, VALIDATION MODE!"

  SYSTEM_REF_PDB="system_ref.pdb"
  log "${REFERENCE_PDB} is copied to: ${SYSTEM_REF_PDB}, will edit it!"
  cp ${REFERENCE_PDB} ${SYSTEM_REF_PDB} #creating local version, we will modify it

  log "Making sure Reference PDB is compatible with mobywat."
  run_step "${SCRIPT_DIR}/pdb-preprocessor.py" --reference ${SYSTEM_REF_PDB}

  # One spec per file: -t selects in the trajectory, the REMARK in the reference.
  MOBYWAT_TARGET="$(mobywat_target_spec "${MOBYWAT_INPUT_PDB}")"
  MOBYWAT_REF_TARGET="$(mobywat_target_spec "${SYSTEM_REF_PDB}")"
  log "MobyWat target: trajectory $MOBYWAT_TARGET, reference $MOBYWAT_REF_TARGET"

  log "Adding mobywat params to ${SYSTEM_REF_PDB}"
  run_step "${SCRIPT_DIR}/add-mobywat-analysis-params.sh" ${SYSTEM_REF_PDB} "$MOBYWAT_REF_TARGET"

  log "Remove TER operator ID if present from ${SYSTEM_REF_PDB}"
  run_step "${SCRIPT_DIR}"/remove-ter-id.sh ${SYSTEM_REF_PDB}

  # Only pass a reference in validation mode.
  MOBYWAT_ARGS+=(-r "${SYSTEM_REF_PDB}")

  log "Running mobywat validation"
else
  log "REFERENCE_PDB is missing or empty, PREDICTION MODE!"
  MOBYWAT_TARGET="$(mobywat_target_spec "${MOBYWAT_INPUT_PDB}")"
  log "MobyWat target: trajectory $MOBYWAT_TARGET (prediction only)"
fi

run_step mobywat "${MOBYWAT_ARGS[@]}" -t "$MOBYWAT_TARGET" -w Auto -n "1-$MOBYWAT_FRAMES" -cls MER -m Prediction -v Diagnostic

if [[ "$MOBYWAT_DEBUG_ENABLED" == true && -n "${REFERENCE_PDB:-}" ]]; then
    MOBYWAT_ANALYSIS_T0=$SECONDS
    if ! run_step mobywat "${MOBYWAT_ARGS[@]}" -t "$MOBYWAT_TARGET" -w Auto -n "1-$MOBYWAT_FRAMES" -m Analysis; then
        log "WARNING: MobyWat Analysis failed; research.sh's sr_frame_* columns will be empty."
    fi
    log "MobyWat Analysis runtime: $(( SECONDS - MOBYWAT_ANALYSIS_T0 )) seconds"
fi

log "tinker.sh completed successfully"