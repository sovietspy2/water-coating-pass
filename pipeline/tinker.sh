#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

INPUT_PDB="$1"
INPUT_ABS="$(cd -P "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"
INPUT_DIR="$(dirname "$INPUT_ABS")"
BASE="$(basename "${INPUT_PDB%.pdb}")"

SECONDS=0
TIMEFILE="${INPUT_ABS}-mm-process-time.txt"

write_runtime() {
  local elapsed="$SECONDS"
  # Write runtime (in seconds) to the requested file
  printf 'runtime for mm is %s seconds\n' "$elapsed" > "$TIMEFILE"
}
trap write_runtime EXIT

pushd "$INPUT_DIR" >/dev/null

pdbxyz "$INPUT_ABS" < "$SCRIPT_DIR/pdbstep.txt"

minimize "$BASE".xyz < "$SCRIPT_DIR/minimizestep.txt"

awk '$3 == "CA" {print "RESTRICT  " $2 "  10"}' "$INPUT_ABS" > key.key

cat "$SCRIPT_DIR/additional-params.key" >> key.key

dynamic "$BASE".xyz_2 -k key.key < "$SCRIPT_DIR/dynamicstep.txt"

xyzpdb "$BASE".arc < "$SCRIPT_DIR/xyzpdbstep.txt"

OUT="${BASE}_mm.pdb"

cp "$BASE".pdb_2 "$OUT"