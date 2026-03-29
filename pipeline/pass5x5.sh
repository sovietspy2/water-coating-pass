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
)

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
      mv -- "$item" "$target"/
    fi
  done
}

cd "$INPUT_DIR"

current_in="$INPUT_PDB"

for i in {1..5}; do
  run_dir="$TOPDIR/$i"
  mkdir -p -- "$run_dir"

  pass_out="${current_in%.pdb}_1WAT.pdb"
  mm_out="${current_in%.pdb}_1WAT_mm.pdb"

  pass "$current_in" 1.8 3.5 1

  if [[ ! -f "$pass_out" ]]; then
    echo "Error: expected pass output not found: $pass_out" >&2
    exit 1
  fi

  # make sure there is no old filtered_renum.pdb present
  rm -f -- filtered_renum.pdb

  # run mm
  run_mm_step "$MODE" "$pass_out" "$SCRIPT_DIR"

  if [[ ! -f "filtered_renum.pdb" ]]; then
    echo "Error: expected MM output not found: filtered_renum.pdb" >&2
    exit 1
  fi

  cp -f -- "filtered_renum.pdb" "$mm_out"

  # KEEP the next input file in the working dir, unless last run
  if [[ $i -ne 5 ]]; then
    WHITELIST+=("$mm_out")
  fi

  move_everything_except_whitelist "$run_dir"

  cp -f -- "$run_dir/$(basename "$mm_out")" .

  current_in="$mm_out"
done

echo "Done. Results are under: $TOPDIR/{1..5}/"