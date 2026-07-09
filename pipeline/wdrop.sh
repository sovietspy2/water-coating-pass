#!/bin/bash
set -euo pipefail

echo "__          _______  _____   ____  _____  "
echo "\ \        / /  __ \|  __ \ / __ \|  __ \ "
echo " \ \  /\  / /| |  | | |__) | |  | | |__) |"
echo "  \ \/  \/ / | |  | |  _  /| |  | |  ___/ "
echo "   \  /\  /  | |__| | | \ \| |__| | |     "
echo "    \/  \/   |_____/|_|  \_\\____/|_|     "

# Usage: ./wdrop.sh INPUT_PDB MODE [REFERENCE_PDB] [--layers L] [--intermediate-md-ns NS] [--final-md-ns NS]
# Example prediciton mode: ./wdrop.sh /abs/path/ASD.pdb gromacs --layers 5
# Example validation mode: ./wdrop.sh /abs/path/ASD.pdb gromacs /abs/path/ASD_crsyst.pdb --layers 5
# Example per-layer MD:    ./wdrop.sh /abs/path/ASD.pdb tinker --intermediate-md-ns 0.1 --final-md-ns 0.5

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
  echo "Usage: $0 <INPUT_PDB> <MODE> [REFERENCE_PDB] [--layers L] [--intermediate-md-ns NS] [--final-md-ns NS]" >&2
  echo "Input PDB: mandatory, used for mobywat prediction mode" >&2
  echo "Modes: gromacs, tinker" >&2
  echo "--layers L: number of water layers (default 5); one layer is deposited+minimized per iteration, so L is also the iteration count" >&2
  echo "--intermediate-md-ns NS: MD length in ns run after each intermediate iteration (default 0); 0 = no intermediate MD (minimize only)" >&2
  echo "--final-md-ns NS: MD length in ns for the final iteration (default 0.1); must be > 0" >&2
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
  local CYCLE_MD_DURATION="$3" # ns for this cycle; >0 runs MD, 0 skips it
  local REFERENCE_PDB_C="$4"
  local RUN_MOBYWAT_C="$5" # 1 on the final iteration (run MobyWat), 0 otherwise

  log "Current input: $INPUT_PDB_LOCAL"
  log "Expected wdrop output: $WDROP_OUTPUT"
  log "Expected mm output: $MM_OUT"

  run_step wdrop --file "${INPUT_PDB_LOCAL}" --sigma 1.8 --weed-dist 3.5 --layers "${WATERS}"
  require_file "$WDROP_OUTPUT" "wdrop output"
  log "wdrop output found: $WDROP_OUTPUT"

  rm -f -- next_step.pdb
  log "Removed old next_step.pdb if present"

  run_step run_mm_step "$MODE" "$WDROP_OUTPUT" "$SCRIPT_DIR" "$CYCLE_MD_DURATION" "$REFERENCE_PDB_C" "$RUN_MOBYWAT_C"
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

# Parse --layers/--intermediate-md-ns/--final-md-ns flags (accepted in any position)
# and the positionals: INPUT_PDB (1), MODE (2), REFERENCE_PDB (3, optional).
LAYERS=5
INTERMEDIATE_MD_NS="0"          # MD length in ns run after each intermediate iteration; 0 = none (minimize only)
FINAL_MD_NS="0.1"               # MD length in ns for the final iteration; must be > 0
POSITIONAL=()
while (( $# )); do
  case "$1" in
    --layers)
      [[ $# -ge 2 ]] || { echo "Error: --layers requires a value" >&2; usage; }
      LAYERS="$2"; shift 2 ;;
    --layers=*)     LAYERS="${1#*=}"; shift ;;
    --intermediate-md-ns)
      [[ $# -ge 2 ]] || { echo "Error: --intermediate-md-ns requires a value" >&2; usage; }
      INTERMEDIATE_MD_NS="$2"; shift 2 ;;
    --intermediate-md-ns=*) INTERMEDIATE_MD_NS="${1#*=}"; shift ;;
    --final-md-ns)
      [[ $# -ge 2 ]] || { echo "Error: --final-md-ns requires a value" >&2; usage; }
      FINAL_MD_NS="$2"; shift 2 ;;
    --final-md-ns=*) FINAL_MD_NS="${1#*=}"; shift ;;
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

# --layers = number of water layers. One layer is deposited + minimized per
# iteration, so LAYERS is also the iteration count (ITERATIONS == LAYERS).
[[ "$LAYERS" =~ ^[1-9][0-9]*$ ]] || {
  echo "Error: --layers must be a positive integer (got '$LAYERS')" >&2; exit 1; }

# --intermediate-md-ns = MD length in ns run after each intermediate iteration.
# Non-negative decimal; 0 = no intermediate MD (deposit + minimize only). This single
# knob decides whether intermediate iterations run MD — there is no separate mode flag.
[[ "$INTERMEDIATE_MD_NS" =~ ^[0-9]+(\.[0-9]+)?$ ]] || {
  echo "Error: --intermediate-md-ns must be a non-negative number in ns (got '$INTERMEDIATE_MD_NS')" >&2; exit 1; }

# --final-md-ns = MD length in ns for the final iteration. Must be a positive
# decimal (e.g. 0.1, 0.001, 1); 0 is rejected because the final iteration always
# runs MD + MobyWat.
{ [[ "$FINAL_MD_NS" =~ ^[0-9]+(\.[0-9]+)?$ ]] && md_enabled "$FINAL_MD_NS"; } || {
  echo "Error: --final-md-ns must be a positive number in ns (got '$FINAL_MD_NS')" >&2; exit 1; }

readonly ITERATIONS="$LAYERS" # one layer deposited per iteration
readonly WATERS_LAYERS_PER_RUN=1
readonly LAYERS INTERMEDIATE_MD_NS FINAL_MD_NS
# Output-dir tag encodes the run parameters: layer count + the two MD lengths (ns).
readonly OUTPUT_TAG="_l${LAYERS}_int${INTERMEDIATE_MD_NS}_fin${FINAL_MD_NS}"

setup_logging "$LOGFILE"

log "Starting wdrop pipeline"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "LAYERS=$LAYERS (one layer per iteration -> $ITERATIONS iterations)"
log "INTERMEDIATE_MD_NS=$INTERMEDIATE_MD_NS ns (intermediate iterations; 0 = none)"
log "FINAL_MD_NS=$FINAL_MD_NS ns (final iteration)"
log "INPUT_DIR=$INPUT_DIR"
log "REFERENCE_PDB(optional)=$REFERENCE_PDB"
validate_script_dir_not_input_dir "$1"

# Step 0B (pdb-preprocessor --target) rewrites the input PDB in place
PDB_NAME="$(basename "$INPUT_PDB" .pdb)"
ORIGINAL_PDB="${INPUT_DIR}/${PDB_NAME}_original.pdb"
log "Step 0A: Backing up untouched input PDB to: $ORIGINAL_PDB"
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

  # MobyWat runs only on the final iteration. The final iteration always runs MD
  # (--final-md-ns); intermediate iterations run MD iff --intermediate-md-ns > 0.
  # The backend gates MD on the duration (>0) and MobyWat on the RUN_MOBYWAT flag.
  if (( I == ITERATIONS )); then
    CYCLE_MD_DURATION=$FINAL_MD_NS
    RUN_MOBYWAT=1
    log "Iteration $I/$ITERATIONS: final iteration — wdrop + minimize + MD ($FINAL_MD_NS ns) + MobyWat"
  else
    CYCLE_MD_DURATION=$INTERMEDIATE_MD_NS
    RUN_MOBYWAT=0
    if md_enabled "$INTERMEDIATE_MD_NS"; then
      log "Iteration $I/$ITERATIONS: intermediate iteration — wdrop + minimize + MD ($INTERMEDIATE_MD_NS ns), no MobyWat"
    else
      log "Iteration $I/$ITERATIONS: intermediate iteration — wdrop + minimize only, no MD"
    fi
  fi

  if (( ITERATIONS > 1 )); then
    RUN_DIR="$OUTPUT_ROOT/$I"
  else
    RUN_DIR="$OUTPUT_ROOT"
  fi
  mkdir -p -- "$RUN_DIR"
  log "Run directory: $RUN_DIR"

  run_wdrop_and_mm "$CURRENT_IN" "$WATERS_LAYERS_PER_RUN" "$CYCLE_MD_DURATION" "$REFERENCE_PDB" "$RUN_MOBYWAT"

  move_everything_except "$RUN_DIR" "${BASE_WHITELIST[@]}"
  require_file "$RUN_DIR/$(basename "$LAST_MM_OUT")" "moved MM result"

  if (( I < ITERATIONS )); then
    NEXT_IN="cycle_input.pdb"
    cp -f -- "$RUN_DIR/$(basename "$LAST_MM_OUT")" "./$NEXT_IN"
    # Needed for tinker output pdb, tinker xyzpdb is sometimes using more columns than PDB standard defines
    run_step "$SCRIPT_DIR/format_pdb.py" "./$NEXT_IN"
    CURRENT_IN="$NEXT_IN"
    log "Next iteration input set to: $CURRENT_IN (columns normalized; copied from $(basename "$LAST_MM_OUT"))"
  fi
done

if (( ITERATIONS > 1 )); then
  log "Done. Results are under: $OUTPUT_ROOT/{1..$ITERATIONS}/"
else
  log "Done. Results are under: $OUTPUT_ROOT"
fi
