#!/bin/bash
set -euo pipefail

# Usage: ./run_pipeline.sh INPUT_PDB MODE
# Example: ./run_pipeline.sh /abs/path/ASD.pdb gromacs

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
WAT_PDB="${INPUT_DIR}/${BASE}_5WAT.pdb"

LOGFILE="${INPUT_DIR}/application.LOG"
: > "$LOGFILE"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGFILE"
}

run_step() {
  log "RUN: $*"
  "$@"
  local rc=$?
  log "DONE (rc=$rc): $*"
  return "$rc"
}

trap 'rc=$?; log "FAILED at line $LINENO with exit code $rc"; exit $rc' ERR

log "Starting run_pipeline"
log "INPUT_PDB=$INPUT_PDB"
log "MODE=$MODE"
log "INPUT_DIR=$INPUT_DIR"
log "WAT_PDB=$WAT_PDB"

# ---- 1) Run pass----
log "Step 1: running pass"
run_step pass "$INPUT_PDB" 1.8 3.5 5

if [[ ! -f "$WAT_PDB" ]]; then
  log "ERROR: expected pass output not found: $WAT_PDB"
  exit 1
fi
log "Pass output found: $WAT_PDB"

# ---- 2) Run minimization once with INPUT_PDB_5WAT ----
log "Step 2: running minimization"
run_step run_mm_step "$MODE" "$WAT_PDB" "$SCRIPT_DIR"

# ---- 3)
if [[ ! -f "${INPUT_DIR}/filtered_renum.pdb" ]]; then
  log "ERROR: filtered_renum.pdb not found after minimization"
  exit 1
fi

WAT_BASE="$(basename "${WAT_PDB%.pdb}")"
log "Creating new file: ${WAT_BASE}-mm.pdb"
cp -- "${INPUT_DIR}/filtered_renum.pdb" "${INPUT_DIR}/${WAT_BASE}-mm.pdb"
log "Copied filtered_renum.pdb -> ${INPUT_DIR}/${WAT_BASE}-mm.pdb"

# ---- 4) Move outputs into a folder named after INPUT_PDB (base) ----
OUTDIR="$(make_output_dir "$INPUT_PDB" "5x1")"
mkdir -p -- "$OUTDIR"
log "Output directory: $OUTDIR"

# Fill this whitelist later (exact filenames or glob patterns).
# These names are excluded from moving.
WHITELIST=(
  "$BASE.pdb"
  "$(basename "$OUTDIR")"
)

log "Whitelist entries: ${WHITELIST[*]}"

# Use extglob to move everything except the whitelist
shopt -s extglob nullglob
pattern="!("
for w in "${WHITELIST[@]}"; do
  pattern+="$w|"
done
pattern="${pattern%|})"

log "Move pattern: $pattern"

# Avoid moving the output directory into itself
cd "$INPUT_DIR"
log "Changed working directory to: $INPUT_DIR"

log "Moving files into: $OUTDIR"
mv -- $pattern "$OUTDIR"/
log "Move completed"

log "Done. Results are under: $OUTDIR"
```
