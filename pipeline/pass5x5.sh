#!/bin/bash
set -euo pipefail

# Usage: ./pass5x5.sh INPUT_PDB MODE
# Example: ./pass5x5.sh /abs/path/ASD.pdb gromacs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <INPUT_PDB> <MODE>" >&2
  echo "Modes: gromacs, tinker" >&2
  exit 1
fi

INPUT_PDB="$(normalize_input_pdb "$1")"
MODE="$2"

if [[ "$MODE" != "gromacs" && "$MODE" != "tinker" ]]; then
  echo "Error: invalid mode '$MODE' (expected gromacs or tinker)" >&2
  exit 1
fi

INPUT_DIR="$(dirname "$INPUT_PDB")"
BASE="$(basename "${INPUT_PDB%.pdb}")"

TOPDIR="$(make_output_dir "$INPUT_PDB" "5x5")"
mkdir -p -- "$TOPDIR"

WHITELIST=(
  "$(basename "$INPUT_PDB")"
  "$(basename "$TOPDIR")"
  "application.LOG"
)

LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging "$LOGFILE"

trap 'rc=$?; log "FAILED at line $LINENO with exit code $rc"; exit $rc' ERR

move_everything_except_whitelist() {
  local target="$1"
  local item
  local skip

  for item in ./*; do
    [[ -e "$item" ]] || continue
    item="${item#./}"
    skip=0

    for w in "${WHITELIST[@]}"; do
      if [[ "$item" == "$w" ]]; then
        skip=1
        break
      fi
    done

    if [[ "$skip" -eq 0 ]]; then
      log "Moving: $item -> $target/"
      mv -- "$item" "$target"/
    fi
  done
}

log "Starting pass5x5"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "TOPDIR=$TOPDIR"

cd "$INPUT_DIR"
log "Working directory: $INPUT_DIR"

log "Step 0: Removing multiple models from PDB, keeping the first one."
run_step "$SCRIPT_DIR"/model-reducer.sh "$INPUT_PDB"

current_in="$INPUT_PDB"

for i in {1..5}; do
  log "===== ITERATION $i/5 ====="
  run_dir="$TOPDIR/$i"
  mkdir -p -- "$run_dir"
  log "Run directory: $run_dir"

  pass_out="${current_in%.pdb}_1WAT.pdb"
  mm_out="${current_in%.pdb}_1WAT_mm.pdb"

  log "Current input: $current_in"
  log "Expected pass output: $pass_out"
  log "Expected mm output: $mm_out"

  run_step pass "$current_in" 1.8 3.5 1

  if [[ ! -f "$pass_out" ]]; then
    log "ERROR: expected pass output not found: $pass_out"
    exit 1
  fi
  log "Pass output found: $pass_out"

  # make sure there is no old filtered_renum.pdb present
  rm -f -- filtered_renum.pdb
  log "Removed old filtered_renum.pdb if present"

  # run mm
  run_step run_mm_step "$MODE" "$pass_out" "$SCRIPT_DIR"

  if [[ ! -f "filtered_renum.pdb" ]]; then
    log "ERROR: expected MM output not found: filtered_renum.pdb"
    exit 1
  fi
  log "MM output found: filtered_renum.pdb"

  cp -f -- "filtered_renum.pdb" "$mm_out"
  log "Copied filtered_renum.pdb -> $mm_out"

  # KEEP the next input file in the working dir, unless last run
  if [[ $i -ne 5 ]]; then
    WHITELIST+=("$mm_out")
    log "Added to whitelist for next iteration: $mm_out"
  else
    log "Last iteration: not adding $mm_out to whitelist"
  fi

  move_everything_except_whitelist "$run_dir"

  if [[ ! -f "$run_dir/$(basename "$mm_out")" ]]; then
    log "ERROR: expected file not found after move: $run_dir/$(basename "$mm_out")"
    exit 1
  fi

  cp -f -- "$run_dir/$(basename "$mm_out")" .
  log "Copied next input back into working dir: $(basename "$mm_out")"

  current_in="$mm_out"
  log "Next iteration input set to: $current_in"
done

log "Done. Results are under: $TOPDIR/{1..5}/"