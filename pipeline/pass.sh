#!/bin/bash
set -euo pipefail

# Usage: ./pass.sh INPUT_PDB MODE RUN_TYPE
# Example: ./pass.sh /abs/path/ASD.pdb gromacs LONG

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"
SECONDS=0 # global timer

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <INPUT_PDB> <MODE> <RUN_TYPE>" >&2
  echo "Modes: gromacs, tinker" >&2
  echo "Run types: LONG, SHORT" >&2
  exit 1
fi

INPUT_PDB="$(normalize_input_pdb "$1")"
MODE="$2"
RUN_TYPE="${3^^}"

if [[ "$MODE" != "gromacs" && "$MODE" != "tinker" ]]; then
  echo "Error: invalid mode '$MODE' (expected gromacs or tinker)" >&2
  exit 1
fi

if [[ "$RUN_TYPE" != "LONG" && "$RUN_TYPE" != "SHORT" ]]; then
  echo "Error: invalid RUN_TYPE '$RUN_TYPE' (expected LONG or SHORT)" >&2
  exit 1
fi

INPUT_DIR="$(dirname "$INPUT_PDB")"
BASE="$(basename "${INPUT_PDB%.pdb}")"
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
      mv -- "$item" "$target/"
    fi
  done
}

log "Starting pass pipeline"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "RUN_TYPE=$RUN_TYPE"
log "INPUT_DIR=$INPUT_DIR"

log "Step 0: removing multiple models from PDB, keeping the first one."
run_step "$SCRIPT_DIR"/model-reducer.sh "$INPUT_PDB"

if [[ "$RUN_TYPE" == "SHORT" ]]; then
  WAT_PDB="${INPUT_DIR}/${BASE}_5WAT.pdb"

  log "SHORT flow selected"
  log "WAT_PDB=$WAT_PDB"

  log "Step 1: running pass"
  run_step pass "$INPUT_PDB" 1.8 3.5 5

  if [[ ! -f "$WAT_PDB" ]]; then
    log "ERROR: expected pass output not found: $WAT_PDB"
    exit 1
  fi
  log "Pass output found: $WAT_PDB"

  log "Step 2: running minimization"
  run_step run_mm_step "$MODE" "$WAT_PDB" "$SCRIPT_DIR"

  if [[ ! -f "${INPUT_DIR}/filtered_renum.pdb" ]]; then
    log "ERROR: filtered_renum.pdb not found after minimization"
    exit 1
  fi

  WAT_BASE="$(basename "${WAT_PDB%.pdb}")"
  log "Creating new file: ${WAT_BASE}-mm.pdb"
  cp -- "${INPUT_DIR}/filtered_renum.pdb" "${INPUT_DIR}/${WAT_BASE}-mm.pdb"
  log "Copied filtered_renum.pdb -> ${INPUT_DIR}/${WAT_BASE}-mm.pdb"

  OUTDIR="$(make_output_dir "$INPUT_PDB" "5x1")"
  mkdir -p -- "$OUTDIR"
  log "Output directory: $OUTDIR"

  WHITELIST=(
    "$BASE.pdb"
    "$(basename "$OUTDIR")"
    "application.LOG"
  )

  shopt -s extglob nullglob
  pattern="!("
  for w in "${WHITELIST[@]}"; do
    pattern+="$w|"
  done
  pattern="${pattern%|})"

  log "Whitelist entries: ${WHITELIST[*]}"
  log "Move pattern: $pattern"

  cd "$INPUT_DIR"
  log "Changed working directory to: $INPUT_DIR"

  log "Moving files into: $OUTDIR"
  mv -- $pattern "$OUTDIR"/
  log "Move completed"
  log "Done. Results are under: $OUTDIR"

else
  TOPDIR="$(make_output_dir "$INPUT_PDB" "5x5")"
  mkdir -p -- "$TOPDIR"

  WHITELIST=(
    "$(basename "$INPUT_PDB")"
    "$(basename "$TOPDIR")"
    "application.LOG"
  )

  log "LONG flow selected"
  log "TOPDIR=$TOPDIR"

  cd "$INPUT_DIR"
  log "Working directory: $INPUT_DIR"

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

    rm -f -- filtered_renum.pdb
    log "Removed old filtered_renum.pdb if present"

    run_step run_mm_step "$MODE" "$pass_out" "$SCRIPT_DIR"

    if [[ ! -f "filtered_renum.pdb" ]]; then
      log "ERROR: expected MM output not found: filtered_renum.pdb"
      exit 1
    fi
    log "MM output found: filtered_renum.pdb"

    cp -f -- "filtered_renum.pdb" "$mm_out"
    log "Copied filtered_renum.pdb -> $mm_out"

    if [[ $i -ne 5 ]]; then
      WHITELIST+=("$(basename "$mm_out")")
      log "Added to whitelist for next iteration: $(basename "$mm_out")"
    else
      log "Last iteration: not adding $(basename "$mm_out") to whitelist"
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
fi

elapsed="$SECONDS"
log "Total runtime was $elapsed seconds"