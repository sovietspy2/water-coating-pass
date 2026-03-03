#!/bin/bash
set -euo pipefail

# Usage: ./pass5x5.sh INPUT_PDB
# Example: ./pass5x5.sh ASD.pdb

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <INPUT_PDB>" >&2
  exit 1
fi

INPUT_PDB="$1"
if [[ ! -f "$INPUT_PDB" ]]; then
  echo "Error: file not found: $INPUT_PDB" >&2
  exit 1
fi

BASE="${INPUT_PDB%.pdb}"
TOPDIR="${BASE}5x5"
mkdir -p -- "$TOPDIR"

WHITELIST=(
  "pass5x1.sh"
  "pass5x5.sh"
  "mm_minim.sh"
  "pass"
  "cg.mdp"
  "md.mdp"
  "clean_pdb.py"
  "st.mdp"
  "cleanup.sh"
  "$INPUT_PDB"
  "$TOPDIR"
)

move_everything_except_whitelist() {
  local target="$1"
  shopt -s extglob nullglob

  local pattern="!("
  for w in "${WHITELIST[@]}"; do
    pattern+="$w|"
  done
  pattern="${pattern%|})"

  mv -- $pattern "$target"/
}

current_in="$INPUT_PDB"

for i in {1..5}; do
  run_dir="$TOPDIR/$i"
  mkdir -p -- "$run_dir"

  # pass output naming rule:
  # ASD.pdb -> ASD_1WAT.pdb, then ASD_1WAT.pdb -> ASD_1WAT_1WAT.pdb, ...
  pass_out="${current_in%.pdb}_1WAT.pdb"

  # mm_minim input naming rule (your examples):
  # ASD.pdb -> ASD_WAT1.pdb
  # ASD_1WAT.pdb -> ASD_WAT1_1WAT.pdb
  # i.e., take pass_out and replace the first occurrence of "_1WAT" with "_WAT1"
  mm_in="${pass_out/_1WAT/_WAT1}"   # bash string replacement: ${var/pattern/repl} [web:81]

  echo "[run $i] pass \"$current_in\" 1.8 3.5 1  (expect: $pass_out)"
  pass "$current_in" 1.8 3.5 1

  if [[ ! -f "$pass_out" ]]; then
    echo "Error: expected pass output not found: $pass_out" >&2
    exit 1
  fi

  # If pass already produces mm_in directly, this will be a no-op overwrite-safe copy.
  # If it does NOT, we create mm_in from pass_out so mm_minim.sh matches your required name.
  if [[ "$mm_in" != "$pass_out" ]]; then
    cp -f -- "$pass_out" "$mm_in"
  fi

  echo "[run $i] mm_minim.sh \"$mm_in\""
  ./mm_minim.sh "$mm_in"

  # Move run outputs into the numbered folder
  move_everything_except_whitelist "$run_dir"

  # Prepare next iteration input (bring it back into working dir)
  cp -f -- "$run_dir/$pass_out" .
  current_in="$pass_out"
done

echo "Done. Results are under: $TOPDIR/{1..5}/"
