#!/bin/bash
set -euo pipefail

echo "__          _______  _____   ____  _____  "
echo "\ \        / /  __ \|  __ \ / __ \|  __ \ "
echo " \ \  /\  / /| |  | | |__) | |  | | |__) |"
echo "  \ \/  \/ / | |  | |  _  /| |  | |  ___/ "
echo "   \  /\  /  | |__| | | \ \| |__| | |     "
echo "    \/  \/   |_____/|_|  \_\\____/|_|     "

# Usage: ./wdrop.sh INPUT_PDB MODE [REFERENCE_PDB] [--iterations N] [--layers L]
# Example prediciton mode: ./wdrop.sh /abs/path/ASD.pdb gromacs --iterations 5
# Example validation mode: ./wdrop.sh /abs/path/ASD.pdb gromacs /abs/path/ASD_crsyst.pdb --iterations 5

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
  echo "Usage: $0 <INPUT_PDB> <MODE> [REFERENCE_PDB] [--iterations N] [--layers L]" >&2
  echo "Input PDB: mandatory, used for mobywat prediction mode" >&2
  echo "Modes: gromacs, tinker" >&2
  echo "--iterations N: number of deposit+minimize cycles (default 1)" >&2
  echo "--layers L: total water layers across the run (default 5); each cycle deposits L/N" >&2
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

# Parse --iterations/--layers flags (accepted in any position) and the
# positionals: INPUT_PDB (1), MODE (2), REFERENCE_PDB (3, optional).
ITERATIONS=1
LAYERS=5
POSITIONAL=()
while (( $# )); do
  case "$1" in
    --iterations)
      [[ $# -ge 2 ]] || { echo "Error: --iterations requires a value" >&2; usage; }
      ITERATIONS="$2"; shift 2 ;;
    --iterations=*) ITERATIONS="${1#*=}"; shift ;;
    --layers)
      [[ $# -ge 2 ]] || { echo "Error: --layers requires a value" >&2; usage; }
      LAYERS="$2"; shift 2 ;;
    --layers=*)     LAYERS="${1#*=}"; shift ;;
    -h|--help)      usage ;;
    --)             shift; while (( $# )); do POSITIONAL+=("$1"); shift; done ;;
    -*)             echo "Error: unknown option '$1'" >&2; usage ;;
    *)              POSITIONAL+=("$1"); shift ;;
  esac
done
if (( ${#POSITIONAL[@]} > 0 )); then
  set -- "${POSITIONAL[@]}"
else
  set --
fi

[[ $# -eq 2 || $# -eq 3 ]] || usage

INPUT_PDB="$(normalize_input_pdb "$1")"
readonly INPUT_PDB
readonly MODE="$2"
readonly REFERENCE_PDB="${3:-}"
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

# --iterations = number of deposit+minimize cycles; --layers = TOTAL water layers
# across the run. Each cycle deposits LAYERS/ITERATIONS layers via the wdrop binary,
# so LAYERS must be an exact multiple of ITERATIONS.
[[ "$ITERATIONS" =~ ^[1-9][0-9]*$ ]] || {
  echo "Error: --iterations must be a positive integer (got '$ITERATIONS')" >&2; exit 1; }
[[ "$LAYERS" =~ ^[1-9][0-9]*$ ]] || {
  echo "Error: --layers must be a positive integer (got '$LAYERS')" >&2; exit 1; }
(( LAYERS % ITERATIONS == 0 )) || {
  echo "Error: --layers ($LAYERS) must be divisible by --iterations ($ITERATIONS)" >&2; exit 1; }

readonly ITERATIONS LAYERS
readonly WATERS_LAYERS_PER_RUN=$(( LAYERS / ITERATIONS ))
readonly OUTPUT_TAG="_i${ITERATIONS}_l${LAYERS}"
# MD + MobyWat run only on the final cycle, so a single MD duration is needed.
readonly FINAL_MD_DURATION="0.001" # in ns

setup_logging "$LOGFILE"

log "Starting wdrop pipeline"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "ITERATIONS=$ITERATIONS"
log "LAYERS=$LAYERS (total; $WATERS_LAYERS_PER_RUN per cycle)"
log "INPUT_DIR=$INPUT_DIR"
log "REFERENCE_PDB(optional)=$REFERENCE_PDB"
validate_script_dir_not_input_dir "$1"

# The script is going to change the input PDB, we are making a reference file from the original input for MobyWat validation purposes
PDB_NAME="$(basename "$INPUT_PDB" .pdb)"
ORIGINAL_PDB="${INPUT_DIR}/${PDB_NAME}_original.pdb"
log "Step 0A: Creating reference PDB file: $ORIGINAL_PDB)"
cp "$INPUT_PDB" "$ORIGINAL_PDB"

# This step is handling X-ray crystallography created PDB-s, adding residues if necessary, also it removes existing waters
log "Step 0B: X-ray PDB fix, fixing missing atoms, removing existing waters, reducing to single model! PYTHON DEPENDENCY!"
run_step "$SCRIPT_DIR/pdb-preprocessor.py" --target "$INPUT_PDB"

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
    log "Iteration $I/$ITERATIONS: final iteration — wdrop + minimize + MD ($FINAL_MD_DURATION ns) + MobyWat"
  else
    MOBYWAT_OUTPUT_ENABLED=false
    log "Iteration $I/$ITERATIONS: intermediate iteration — wdrop + minimize only, no MD"
  fi
  # MD only runs on the final iteration (gated by MOBYWAT_OUTPUT_ENABLED in the backend);
  # MD_DURATION is ignored when MobyWat output is disabled.
  MD_DURATION=$FINAL_MD_DURATION

  if (( ITERATIONS > 1 )); then
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

if (( ITERATIONS > 1 )); then
  log "Done. Results are under: $OUTPUT_ROOT/{1..$ITERATIONS}/"
else
  log "Done. Results are under: $OUTPUT_ROOT"
fi