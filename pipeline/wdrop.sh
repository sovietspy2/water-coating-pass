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

# Python VENV setup
VENV_DIR=".venv"

if [ ! -d "$VENV_DIR" ]; then
    echo "[INFO] Creating virtual environment..."
    python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

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
  local path="$1"
  local label="$2"

  if [[ ! -f "$path" ]]; then
    log "ERROR: expected $label not found: $path"
    exit 1
  fi
}

move_everything_except() {
  local target="$1"
  shift
  local whitelist=("$@")
  local item
  local keep
  local name

  mkdir -p -- "$target"

  for item in ./*; do
    [[ -e "$item" ]] || continue
    name="${item#./}"
    keep=0

    for w in "${whitelist[@]}"; do
      if [[ "$name" == "$w" ]]; then
        keep=1
        break
      fi
    done

    if [[ "$keep" -eq 0 ]]; then
      log "Moving: $name -> $target/"
      mv -- "$name" "$target/"
    fi
  done
}

run_wdrop_and_mm() {
  local input_pdb="$1"
  local waters="$2"
  local wdrop_output="${input_pdb%.pdb}_${waters}WAT.pdb"
  local mm_out="${wdrop_output%.pdb}_mm.pdb"
  local md_duration="$3"
  local MOBYWAT_OUTPUT_ENABLED="$4"
  local REFERENCE_PDB_C="$5"

  log "Current input: $input_pdb"
  log "Expected wdrop output: $wdrop_output"
  log "Expected mm output: $mm_out"

  run_step wdrop "$input_pdb" 1.8 3.5 "$waters"
  require_file "$wdrop_output" "wdrop output"
  log "wdrop output found: $wdrop_output"

  rm -f -- next_step.pdb
  log "Removed old next_step.pdb if present"

  run_step run_mm_step "$MODE" "$wdrop_output" "$SCRIPT_DIR" "$md_duration" "$MOBYWAT_OUTPUT_ENABLED" "$REFERENCE_PDB_C"
  require_file "next_step.pdb" "MM output"
  log "MM output found: next_step.pdb"

  cp -f -- next_step.pdb "$mm_out"
  log "Copied next_step.pdb -> $mm_out"

  LAST_MM_OUT="$mm_out"
}

on_error() {
  local rc="$?"
  local line="${BASH_LINENO[0]:-$LINENO}"
  log "FAILED at line $line with exit code $rc"
  exit "$rc"
}

on_exit() {
  log "Total runtime was $SECONDS seconds"
}

trap on_error ERR
trap on_exit EXIT

[[ $# -eq 3 || $# -eq 4 ]] || usage

readonly INPUT_PDB="$(normalize_input_pdb "$1")"
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
run_step "$SCRIPT_DIR/pdb-atom-fixes.py" "$INPUT_PDB"

cd "$INPUT_DIR"
log "Working directory: $INPUT_DIR"

readonly OUTPUT_ROOT="$(make_output_dir "$INPUT_PDB" "$OUTPUT_TAG")"
mkdir -p -- "$OUTPUT_ROOT"
log "Output directory: $OUTPUT_ROOT"

base_whitelist=(
  "$INPUT_FILE"
  "$(basename "$OUTPUT_ROOT")"
  "application.LOG"
  "test.log"
  "result.log"
  "$(basename "$ORIGINAL_PDB")"
  "$(basename "$REFERENCE_PDB")" #Reference PDB might be in the working dir
)

current_in="$INPUT_FILE"

for ((i = 1; i <= ITERATIONS; i++)); do
  log "===== ITERATION $i/$ITERATIONS ====="

  if (( i == ITERATIONS )); then
    MOBYWAT_OUTPUT_ENABLED=true
    MD_DURATION=$FINAL_MD_DURATION
    log "Iteration $i/$ITERATIONS: last iteration, md will run for $MD_DURATION nanosec"
  else
    MOBYWAT_OUTPUT_ENABLED=false
    MD_DURATION=$INTERMEDIATE_MD_DURATION
    log "Iteration $i/$ITERATIONS: not last iteration, md will run for $MD_DURATION nanosec"
  fi

  if [[ "$RUN_TYPE" == "LONG" ]]; then
    run_dir="$OUTPUT_ROOT/$i"
  else
    run_dir="$OUTPUT_ROOT"
  fi
  mkdir -p -- "$run_dir"
  log "Run directory: $run_dir"

  run_wdrop_and_mm "$current_in" "$WATERS_LAYERS_PER_RUN" "$MD_DURATION" "$MOBYWAT_OUTPUT_ENABLED" "$REFERENCE_PDB"

  move_everything_except "$run_dir" "${base_whitelist[@]}"
  require_file "$run_dir/$(basename "$LAST_MM_OUT")" "moved MM result"

  if (( i < ITERATIONS )); then
    cp -f -- "$run_dir/$(basename "$LAST_MM_OUT")" .
    log "Copied next input back into working dir: $(basename "$LAST_MM_OUT")"
    current_in="$(basename "$LAST_MM_OUT")"
    log "Next iteration input set to: $current_in"
  fi
done

if [[ "$RUN_TYPE" == "LONG" ]]; then
  log "Done. Results are under: $OUTPUT_ROOT/{1..5}/"
else
  log "Done. Results are under: $OUTPUT_ROOT"
fi