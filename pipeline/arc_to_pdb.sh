#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "$SCRIPT_DIR/pipeline_common.sh" ]]; then
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/pipeline_common.sh"
fi

if ! declare -F log >/dev/null 2>&1; then
  log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
  }
fi

if ! declare -F run_step >/dev/null 2>&1; then
  run_step() {
    "$@"
  }
fi

usage() {
  cat >&2 <<'EOF'
Usage:
  arc_to_multimodel_pdb.sh <trajectory.arc> [output.pdb] [parameter.prm]

Arguments:
  trajectory.arc   TINKER trajectory archive
  output.pdb       Optional output path (default: <input_dir>/<base>-trajectory.pdb)
  parameter.prm    Optional parameter file (default: <input_dir>/amber99.prm)

Environment:
  KEEP_TMP=1       Keep temporary split-frame files for debugging
EOF
  exit 1
}

[[ $# -ge 1 ]] || usage

ARC="$1"
[[ -f "$ARC" ]] || { echo "ERROR: ARC file not found: $ARC" >&2; exit 1; }

INPUT_DIR="$(cd "$(dirname "$ARC")" && pwd)"
ARC_FILE="$(basename "$ARC")"
PDB_NAME="${ARC_FILE%.arc}"

OUT="${2:-${INPUT_DIR}/${PDB_NAME}-trajectory.pdb}"
PRM="${3:-${INPUT_DIR}/amber99.prm}"
SEQ="${INPUT_DIR}/${PDB_NAME}.seq"

[[ -f "$PRM" ]] || { echo "ERROR: Parameter file not found: $PRM" >&2; exit 1; }
[[ -f "$SEQ" ]] || {
  echo "ERROR: Sequence file not found: $SEQ" >&2
  echo "xyzpdb needs the .seq produced during the original pdbxyz step." >&2
  exit 1
}

LOGFILE="${INPUT_DIR}/application.LOG"
if declare -F setup_logging >/dev/null 2>&1; then
  setup_logging "$LOGFILE"
fi

log "Starting arc_to_multimodel_pdb.sh"
log "ARC=$ARC"
log "INPUT_DIR=$INPUT_DIR"
log "PDB_NAME=$PDB_NAME"
log "OUT=$OUT"
log "PRM=$PRM"
log "SEQ=$SEQ"

TMPDIR="$(mktemp -d "${INPUT_DIR}/${PDB_NAME}.arc2pdb.XXXXXX")"

cleanup() {
  if [[ "${KEEP_TMP:-0}" == "1" ]]; then
    log "KEEP_TMP=1, preserving temporary directory: $TMPDIR"
  else
    rm -rf "$TMPDIR"
  fi
}
trap cleanup EXIT

NATOMS="$(awk 'NR==1 {print $1; exit}' "$ARC")"
SECOND_FIRST="$(awk 'NR==2 {print $1; exit}' "$ARC")"
TOTAL_LINES="$(wc -l < "$ARC" | tr -d ' ')"

if [[ -z "${NATOMS:-}" || ! "$NATOMS" =~ ^[0-9]+$ ]]; then
  log "ERROR: Could not determine atom count from ARC: $ARC"
  exit 1
fi

case "$SECOND_FIRST" in
  *[!0-9]* ) LINES_PER_FRAME=$((NATOMS + 2)) ;;
  * )        LINES_PER_FRAME=$((NATOMS + 1)) ;;
esac

if (( TOTAL_LINES % LINES_PER_FRAME != 0 )); then
  log "ERROR: ARC line count ($TOTAL_LINES) is not divisible by lines/frame ($LINES_PER_FRAME)"
  exit 1
fi

NFRAMES=$((TOTAL_LINES / LINES_PER_FRAME))

if (( NFRAMES < 1 )); then
  log "ERROR: No frames detected in ARC: $ARC"
  exit 1
fi

log "Detected NATOMS=$NATOMS"
log "Detected LINES_PER_FRAME=$LINES_PER_FRAME"
log "Detected NFRAMES=$NFRAMES"

log "Splitting ARC into per-frame XYZ files"
awk -v lpf="$LINES_PER_FRAME" -v prefix="$TMPDIR/${PDB_NAME}_frame_" '
{
  frame = int((NR - 1) / lpf) + 1
  file = sprintf("%s%06d.xyz", prefix, frame)
  if (prev != "" && file != prev) close(prev)
  print >> file
  prev = file
}
END {
  if (prev != "") close(prev)
}
' "$ARC"

: > "$OUT"
model=1
shopt -s nullglob

for ((frame=1; frame<=NFRAMES; frame++)); do
  frame_xyz="${TMPDIR}/${PDB_NAME}_frame_$(printf '%06d' "$frame").xyz"
  frame_base="${frame_xyz%.xyz}"

  if [[ ! -f "$frame_xyz" ]]; then
    log "ERROR: Missing split frame XYZ: $frame_xyz"
    exit 1
  fi

  ln -sf "$SEQ" "${frame_base}.seq"

  log "Converting frame $frame/$NFRAMES"
  run_step xyzpdb "$frame_xyz" <<EOF
$PRM
PDB
EOF

  pdb_candidates=( "${frame_base}.pdb"* )
  if [[ ${#pdb_candidates[@]} -eq 0 ]]; then
    log "ERROR: xyzpdb did not produce a PDB for frame file: $frame_xyz"
    exit 1
  fi

  frame_pdb="${pdb_candidates[0]}"

  if (( model == 1 )); then
    awk '
      /^(HEADER|TITLE |COMPND|SOURCE|KEYWDS|EXPDTA|AUTHOR|REMARK|CRYST1)/ {print}
    ' "$frame_pdb" >> "$OUT"
  fi

  printf 'MODEL     %4d\n' "$model" >> "$OUT"

  awk '
    /^(ATOM  |HETATM|ANISOU|TER   )/ {print}
  ' "$frame_pdb" >> "$OUT"

  printf 'ENDMDL\n' >> "$OUT"
  ((model++))
done

echo "END" >> "$OUT"

if [[ ! -s "$OUT" ]]; then
  log "ERROR: Output file is empty: $OUT"
  exit 1
fi

log "Wrote multi-model PDB: $OUT"
log "arc_to_multimodel_pdb.sh completed successfully"