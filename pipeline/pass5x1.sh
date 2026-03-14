#!/bin/bash
set -euo pipefail

# Usage: ./run_pipeline.sh INPUT_PDB
# Example: ./run_pipeline.sh ASD.pdb

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <INPUT_PDB>" >&2
  exit 1
fi

INPUT_PDB="$1"

if [[ ! -f "$INPUT_PDB" ]]; then
  echo "Error: file not found: $INPUT_PDB" >&2
  exit 1
fi

# Derive names
BASE="${INPUT_PDB%.pdb}"     # strips trailing .pdb if present
WAT_PDB="${BASE}_5WAT.pdb"

# ---- 1) Run pass: creates/uses *_5WAT.pdb as you specified ----
# You said: pass INPUT_PDB_5WAT.pdb 1.8 3.5 5
# If pass actually needs the original input too, tell me and I’ll adjust.
pass "$INPUT_PDB" 1.8 3.5 5

# ---- 2) Run minimization once with INPUT_PDB_5WAT ----
./mm_minim.sh "$WAT_PDB"

# ---- 3) Copy filtered_renum.pdb to BASE-FINAL.pdb ----
if [[ ! -f "filtered_renum.pdb" ]]; then
  echo "Error: filtered_renum.pdb not found after minimization" >&2
  exit 1
fi
cp -- "filtered_renum.pdb" "${WAT_PDB}-mm.pdb"

# ---- 4) Move outputs into a folder named after INPUT_PDB (base) ----
OUTDIR="${BASE}5x1"
mkdir -p -- "$OUTDIR"

# Fill this whitelist later (exact filenames or glob patterns).
# These names are excluded from moving.
WHITELIST=(
  "cleanup.sh"
  "pass5x5.sh"
  "pass5x1.sh"
  "mm_minim.sh"
  "cg.mdp"
  "md.mdp"
  "clean_pdb.py"
  "st.mdp"
  "$OUTDIR"
  "$INPUT_PDB"
)

# Use extglob to move everything except the whitelist
shopt -s extglob nullglob
pattern="!("
for w in "${WHITELIST[@]}"; do
  pattern+="$w|"
done
pattern="${pattern%|})"

# Avoid moving the output directory into itself
mv -- $pattern "$OUTDIR"/
