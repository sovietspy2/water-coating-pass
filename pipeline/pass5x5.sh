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

if [[ "$INPUT_PDB" != *.pdb ]]; then
  echo "Error: input must end with .pdb (got: $INPUT_PDB)" >&2
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

  # Safe for spaces/special characters using NUL delimiters. [web:83]
  local -a find_args
  find_args=( -maxdepth 1 -mindepth 1 )

  for w in "${WHITELIST[@]}"; do
    find_args+=( -not -name "$w" )
  done

  find . "${find_args[@]}" -print0 | xargs -0 -I{} mv -- "{}" "$target"/
}

current_in="$INPUT_PDB"

for i in {1..5}; do
  run_dir="$TOPDIR/$i"
  mkdir -p -- "$run_dir"

  # Append _1WAT before the final .pdb every iteration. [web:49][web:46]
  pass_out="${current_in%.pdb}_1WAT.pdb"

  echo "[run $i] pass \"$current_in\" 1.8 3.5 1  (expect: $pass_out)"
  pass "$current_in" 1.8 3.5 1

  if [[ ! -f "$pass_out" ]]; then
    echo "Error: expected pass output not found: $pass_out" >&2
    exit 1
  fi

  echo "[run $i] mm_minim.sh \"$pass_out\""
  ./mm_minim.sh "$pass_out"

  move_everything_except_whitelist "$run_dir"

  cp -f -- "$run_dir/$pass_out" .
  current_in="$pass_out"
done

echo "Done. Results are under: $TOPDIR/{1..5}/"
