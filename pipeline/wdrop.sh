#!/bin/bash
set -euo pipefail

echo "__          _______  _____   ____  _____  "
echo "\ \        / /  __ \|  __ \ / __ \|  __ \ "
echo " \ \  /\  / /| |  | | |__) | |  | | |__) |"
echo "  \ \/  \/ / | |  | |  _  /| |  | |  ___/ "
echo "   \  /\  /  | |__| | | \ \| |__| | |     "
echo "    \/  \/   |_____/|_|  \_\\____/|_|     "

# Usage: ./wdrop.sh INPUT_PDB MODE RUN_TYPE REFERENCE_PDB
# Example prediciton mode: ./wdrop.sh /abs/path/ASD.pdb gromacs LONG
# Example validation mode: ./wdrop.sh /abs/path/ASD.pdb gromacs LONG /abs/path/ASD_crsyst.pdb

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

# Python VENV setup
VENV_DIR="${SCRIPT_DIR}/../.venv"

# Activating venv, if we dont have one it creates one
[ -d "$VENV_DIR" ] || python3 -m venv "$VENV_DIR"
[ -f "$VENV_DIR/bin/activate" ] || {
    echo "[ERROR] Missing activate script: $VENV_DIR/bin/activate" >&2
    exit 1
}

source "$VENV_DIR/bin/activate"

SECONDS=0

usage() {
  echo "Usage: $0 <INPUT_PDB> <MODE> <RUN_TYPE> <REFERENCE_PDB_WITH_WATER>" >&2
  echo "Input PDB: mandatory used for mobywat prediction mode"
  echo "Modes: gromacs, tinker" >&2
  echo "Run types: LONG, SHORT" >&2
  echo "Reference PDB: (optional) for mobywat validaiton mode" >&2
  exit 1
}

require_file() {
  local PATH_VALUE="$1"
  local LABEL="$2"

  if [[ ! -f "$PATH_VALUE" ]]; then
    log "ERROR: expected $LABEL not found: $PATH_VALUE"
    exit 1
  fi
}

move_everything_except() {
  local TARGET="$1"
  shift
  local WHITELIST=("$@")
  local ITEM
  local KEEP
  local NAME

  mkdir -p -- "$TARGET"

  for ITEM in ./*; do
    [[ -e "$ITEM" ]] || continue
    NAME="${ITEM#./}"
    KEEP=0

    for W in "${WHITELIST[@]}"; do
      if [[ "$NAME" == "$W" ]]; then
        KEEP=1
        break
      fi
    done

    if [[ "$KEEP" -eq 0 ]]; then
      log "Moving: $NAME -> $TARGET/"
      mv -- "$NAME" "$TARGET/"
    fi
  done
}

run_wdrop_and_mm() {
  local INPUT_PDB_LOCAL="$1"
  local WATERS="$2"
  local WDROP_OUTPUT="${INPUT_PDB_LOCAL%.pdb}_${WATERS}WAT.pdb"
  local MM_OUT="${WDROP_OUTPUT%.pdb}_mm.pdb"
  local MD_DURATION="$3"
  local MOBYWAT_OUTPUT_ENABLED="$4"
  local REFERENCE_PDB_C="$5"

  log "Current input: $INPUT_PDB_LOCAL"
  log "Expected wdrop output: $WDROP_OUTPUT"
  log "Expected mm output: $MM_OUT"

  run_step wdrop --file "${INPUT_PDB_LOCAL}" --sigma 1.8 --weed-dist 3.5 --layers "${WATERS}"
  require_file "$WDROP_OUTPUT" "wdrop output"
  log "wdrop output found: $WDROP_OUTPUT"

  rm -f -- next_step.pdb
  log "Removed old next_step.pdb if present"

  run_step run_mm_step "$MODE" "$WDROP_OUTPUT" "$SCRIPT_DIR" "$MD_DURATION" "$MOBYWAT_OUTPUT_ENABLED" "$REFERENCE_PDB_C"
  require_file "next_step.pdb" "MM output"
  log "MM output found: next_step.pdb"

  cp -f -- next_step.pdb "$MM_OUT"
  log "Copied next_step.pdb -> $MM_OUT"

  LAST_MM_OUT="$MM_OUT"
}

on_error() {
  local RC="$?"
  local LINE="${BASH_LINENO[0]:-$LINENO}"
  log "FAILED at line $LINE with exit code $RC"
  exit "$RC"
}

on_exit() {
  log "Total runtime was $SECONDS seconds"
}

trap on_error ERR
trap on_exit EXIT

[[ $# -eq 3 || $# -eq 4 ]] || usage

INPUT_PDB="$(normalize_input_pdb "$1")"
readonly INPUT_PDB
readonly MODE="$2"
readonly RUN_TYPE="${3^^}"
readonly REFERENCE_PDB="${4:-}"
readonly INPUT_DIR="$(dirname "$INPUT_PDB")"
readonly INPUT_FILE="$(basename "$INPUT_PDB")"
LOGFILE="${INPUT_DIR}/application.LOG"

case "$MODE" in
  gromacs|tinker) ;;
  *)
    echo "Error: invalid mode '$MODE' (expected gromacs or tinker)" >&2
    exit 1
    ;;
esac

case "$RUN_TYPE" in
  LONG)
    readonly ITERATIONS=5
    readonly WATERS_LAYERS_PER_RUN=1
    readonly OUTPUT_TAG="5x5"
    readonly FINAL_MD_DURATION="0.1" #in ns
    readonly INTERMEDIATE_MD_DURATION="0.01" # in ns
    ;;
  SHORT)
    readonly ITERATIONS=1
    readonly WATERS_LAYERS_PER_RUN=5 # used to be 5
    readonly OUTPUT_TAG="5x1"
    readonly FINAL_MD_DURATION="0.1" #in ns
    readonly INTERMEDIATE_MD_DURATION="0.1" # in ns
    ;;
  *)
    echo "Error: invalid RUN_TYPE '$RUN_TYPE' (expected LONG or SHORT)" >&2
    exit 1
    ;;
esac

setup_logging "$LOGFILE"

log "Starting wdrop pipeline"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "RUN_TYPE=$RUN_TYPE"
log "INPUT_DIR=$INPUT_DIR"
log "REFERENCE_PDB(optional)=$REFERENCE_PDB"
validate_script_dir_not_input_dir "$1"

# The script is going to change the input PDB, we are making a reference file from the original input for MobyWat validation purposes
PDB_NAME="$(basename "$INPUT_PDB" .pdb)"
ORIGINAL_PDB="${INPUT_DIR}/${PDB_NAME}_original.pdb"
log "Step 0A: Creating reference PDB file: $ORIGINAL_PDB)"
cp "$INPUT_PDB" "$ORIGINAL_PDB"

log "Step 0B: removing multiple models from PDB, keeping the first one."
run_step "$SCRIPT_DIR"/model-reducer.sh "$INPUT_PDB"

# This step is handling X-ray crystallography created PDB-s, adding residues if necessary, also it removes existing waters
log "Step 0C: X-ray PDB fix, fixing missing atoms, removing existing waters! PYTHON DEPENDENCY!"
run_step "$SCRIPT_DIR/target-pdb-preprocessor.py" "$INPUT_PDB"

cd "$INPUT_DIR"
log "Working directory: $INPUT_DIR"

OUTPUT_ROOT="$(make_output_dir "$INPUT_PDB" "$OUTPUT_TAG")"
readonly OUTPUT_ROOT
mkdir -p -- "$OUTPUT_ROOT"
log "Output directory: $OUTPUT_ROOT"

BASE_WHITELIST=(
  "$INPUT_FILE"
  "$(basename "$OUTPUT_ROOT")"
  "application.LOG"
  "test.log"
  "result.log"
  "$(basename "$ORIGINAL_PDB")"
  "$(basename "$REFERENCE_PDB")" #Reference PDB might be in the working dir
)

CURRENT_IN="$INPUT_FILE"

for ((I = 1; I <= ITERATIONS; I++)); do
  log "===== ITERATION $I/$ITERATIONS ====="

  if (( I == ITERATIONS )); then
    MOBYWAT_OUTPUT_ENABLED=true
    MD_DURATION=$FINAL_MD_DURATION
    log "Iteration $I/$ITERATIONS: last iteration, md will run for $MD_DURATION nanosec"
  else
    MOBYWAT_OUTPUT_ENABLED=false
    MD_DURATION=$INTERMEDIATE_MD_DURATION
    log "Iteration $I/$ITERATIONS: not last iteration, md will run for $MD_DURATION nanosec"
  fi

  if [[ "$RUN_TYPE" == "LONG" ]]; then
    RUN_DIR="$OUTPUT_ROOT/$I"
  else
    RUN_DIR="$OUTPUT_ROOT"
  fi
  mkdir -p -- "$RUN_DIR"
  log "Run directory: $RUN_DIR"

  run_wdrop_and_mm "$CURRENT_IN" "$WATERS_LAYERS_PER_RUN" "$MD_DURATION" "$MOBYWAT_OUTPUT_ENABLED" "$REFERENCE_PDB"

  move_everything_except "$RUN_DIR" "${BASE_WHITELIST[@]}"
  require_file "$RUN_DIR/$(basename "$LAST_MM_OUT")" "moved MM result"

  if (( I < ITERATIONS )); then
    cp -f -- "$RUN_DIR/$(basename "$LAST_MM_OUT")" .
    log "Copied next input back into working dir: $(basename "$LAST_MM_OUT")"
    CURRENT_IN="$(basename "$LAST_MM_OUT")"
    log "Next iteration input set to: $CURRENT_IN"
  fi
done

if [[ "$RUN_TYPE" == "LONG" ]]; then
  log "Done. Results are under: $OUTPUT_ROOT/{1..5}/"
else
  log "Done. Results are under: $OUTPUT_ROOT"
fi