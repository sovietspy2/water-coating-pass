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

# ---- 1) Run pass----
pass "$INPUT_PDB" 1.8 3.5 5

# ---- 2) Run minimization once with INPUT_PDB_5WAT ----
run_mm_step "$MODE" "$WAT_PDB" "$SCRIPT_DIR"

# ---- 3)
if [[ ! -f "${INPUT_DIR}/filtered_renum.pdb" ]]; then
  echo "Error: filtered_renum.pdb not found after minimization" >&2
  exit 1
fi

WAT_BASE="$(basename "${WAT_PDB%.pdb}")"
echo "Creating new file: ${WAT_BASE}-mm.pdb"
cp -- "${INPUT_DIR}/filtered_renum.pdb" "${INPUT_DIR}/${WAT_BASE}-mm.pdb"

# ---- 4) Move outputs into a folder named after INPUT_PDB (base) ----
OUTDIR="$(make_output_dir "$INPUT_PDB" "5x1")"
mkdir -p -- "$OUTDIR"

# Fill this whitelist later (exact filenames or glob patterns).
# These names are excluded from moving.
WHITELIST=(
  "pass5x1.sh"
  "pass5x5.sh"
  "pipeline_common.sh"
  # gromacs specific
  "cleanup.sh"
  "mm_minim.sh"
  "gromacs-cg.mdp"
  "gromacs-md.mdp"
  "clean_pdb.py"
  "gromacs-st.mdp"
  # tinker specific
  "amber99.prm"
  "tinker.sh"
  "$BASE.pdb"
  "$(basename "$OUTDIR")"
)

# Use extglob to move everything except the whitelist
shopt -s extglob nullglob
pattern="!("
for w in "${WHITELIST[@]}"; do
  pattern+="$w|"
done
pattern="${pattern%|})"

# Avoid moving the output directory into itself
cd "$INPUT_DIR"
mv -- $pattern "$OUTDIR"/