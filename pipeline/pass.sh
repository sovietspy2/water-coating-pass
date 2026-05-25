#!/bin/bash
set -euo pipefail

# Usage: ./pass.sh INPUT_PDB MODE RUN_TYPE
# Example: ./pass.sh /abs/path/ASD.pdb gromacs LONG

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

SECONDS=0

usage() {
  echo "Usage: $0 <INPUT_PDB> <MODE> <RUN_TYPE>" >&2
  echo "Modes: gromacs, tinker" >&2
  echo "Run types: LONG, SHORT" >&2
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

run_pass_and_mm() {
  local input_pdb="$1"
  local waters="$2"
  local pass_out="${input_pdb%.pdb}_${waters}WAT.pdb"
  local mm_out="${pass_out%.pdb}_mm.pdb"

  log "Current input: $input_pdb"
  log "Expected pass output: $pass_out"
  log "Expected mm output: $mm_out"

  run_step pass "$input_pdb" 1.8 3.5 "$waters"
  require_file "$pass_out" "pass output"
  log "Pass output found: $pass_out"

  rm -f -- filtered_renum.pdb
  log "Removed old filtered_renum.pdb if present"

  run_step run_mm_step "$MODE" "$pass_out" "$SCRIPT_DIR"
  require_file "filtered_renum.pdb" "MM output"
  log "MM output found: filtered_renum.pdb"

  cp -f -- filtered_renum.pdb "$mm_out"
  log "Copied filtered_renum.pdb -> $mm_out"

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

[[ $# -eq 3 ]] || usage

readonly INPUT_PDB="$(normalize_input_pdb "$1")"
readonly MODE="$2"
readonly RUN_TYPE="${3^^}"
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
    readonly WATERS_PER_PASS=1
    readonly OUTPUT_TAG="5x5"
    readonly FINAL_MD_DURATION=1 #in ns
    readonly INTERMEDIATE_MD_DURATION="0.1" # in ns
    ;;
  SHORT)
    readonly ITERATIONS=1
    readonly WATERS_PER_PASS=5
    readonly OUTPUT_TAG="5x1"
    readonly FINAL_MD_DURATION=1 #in ns
    readonly INTERMEDIATE_MD_DURATION=1 # in ns
    ;;
  *)
    echo "Error: invalid RUN_TYPE '$RUN_TYPE' (expected LONG or SHORT)" >&2
    exit 1
    ;;
esac

setup_logging "$LOGFILE"

log "Starting pass pipeline"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "RUN_TYPE=$RUN_TYPE"
log "INPUT_DIR=$INPUT_DIR"
validate_script_dir_not_input_dir "$1"

log "Step 0: removing multiple models from PDB, keeping the first one."
run_step "$SCRIPT_DIR"/model-reducer.sh "$INPUT_PDB"

cd "$INPUT_DIR"
log "Working directory: $INPUT_DIR"

readonly OUTPUT_ROOT="$(make_output_dir "$INPUT_PDB" "$OUTPUT_TAG")"
mkdir -p -- "$OUTPUT_ROOT"
log "Output directory: $OUTPUT_ROOT"

base_whitelist=(
  "$INPUT_FILE"
  "$(basename "$OUTPUT_ROOT")"
  "application.LOG"
)

current_in="$INPUT_FILE"

for ((i = 1; i <= ITERATIONS; i++)); do
  log "===== ITERATION $i/$ITERATIONS ====="

  if (( i == ITERATIONS )); then
    MD_DURATION=$FINAL_MD_DURATION
    log "Iteration $i/$ITERATIONS: last iteration, md will run for $MD_DURATION nanosec"
  else
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

  run_pass_and_mm "$current_in" "$WATERS_PER_PASS" "$MD_DURATION"

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