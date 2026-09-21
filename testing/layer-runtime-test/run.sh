#!/usr/bin/env bash
set -euo pipefail

# Per-layer wdrop runtime benchmark.
# For every PDB in INPUT_DIR: run wdrop with --layers 1, feed its output back as the
# next input, repeat LAYERS times, and record the runtime (ms) of each run as reported
# by wdrop's own log file ("runtime: X ms"). One CSV row per PDB.

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

WDROP_BIN="$PROJECT_ROOT/src/wdrop"
PREPROCESSOR="$PROJECT_ROOT/pipeline/pdb-preprocessor.py"
LAYERS=5
SIGMA_P="1.8"
WEED_DIST="3.5"
OUTPUT_CSV=""
WORK_DIR=""
KEEP_WORK=false
PREPROCESS=false
INPUT_DIR=""

usage() {
  cat <<EOF
Usage: $(basename "$0") <INPUT_DIR> [OPTIONS]

Runs wdrop --layers 1 on every *.pdb in INPUT_DIR, then re-runs it on its own output,
LAYERS times in total, and writes the runtime of each run (ms, from wdrop's log) to a CSV:

  pdb_name,layer1_ms,layer2_ms,layer3_ms,layer4_ms,layer5_ms

The input files are never modified (each PDB is copied into a work directory).
A PDB whose run fails gets FAILED in that column and empty columns after it.

Options:
  -n, --layers N       number of chained runs (default: 5)
  -o, --output FILE    CSV path (default: results/<UTC-timestamp>.csv next to this script)
  -p, --wdrop-bin PATH wdrop executable (default: src/wdrop)
  --preprocess         run pipeline/pdb-preprocessor.py --target on each copy first
                       (removes waters/heterogens, collapses altlocs, rebuilds missing
                       atoms), i.e. the same input wdrop gets in the pipeline's 1st cycle
  --sigma VALUE        probe radius (default: 1.8)
  --weed-dist VALUE    minimal O-O distance (default: 3.5)
  -w, --work-dir DIR   work directory (default: mktemp)
  --keep               keep the work directory (default: removed on exit)
  -h, --help           show this help

Note: wdrop rejects PDBs containing alternate locations; use --preprocess for raw RCSB files.
EOF
}

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -n|--layers)    [[ $# -ge 2 ]] || die "$1 requires a value"; LAYERS="$2"; shift 2 ;;
      -o|--output)    [[ $# -ge 2 ]] || die "$1 requires a value"; OUTPUT_CSV="$2"; shift 2 ;;
      -p|--wdrop-bin) [[ $# -ge 2 ]] || die "$1 requires a value"; WDROP_BIN="$2"; shift 2 ;;
      -w|--work-dir)  [[ $# -ge 2 ]] || die "$1 requires a value"; WORK_DIR="$2"; shift 2 ;;
      --sigma)        [[ $# -ge 2 ]] || die "$1 requires a value"; SIGMA_P="$2"; shift 2 ;;
      --weed-dist)    [[ $# -ge 2 ]] || die "$1 requires a value"; WEED_DIST="$2"; shift 2 ;;
      --preprocess)   PREPROCESS=true; shift ;;
      --keep)         KEEP_WORK=true; shift ;;
      -h|--help)      usage; exit 0 ;;
      -*)             die "unknown option: $1" ;;
      *)
        [[ -z "$INPUT_DIR" ]] || die "only one INPUT_DIR is accepted (got '$INPUT_DIR' and '$1')"
        INPUT_DIR="$1"; shift ;;
    esac
  done
}

cleanup() {
  if [[ -n "${WORK_DIR:-}" && -d "$WORK_DIR" ]]; then
    if [[ "$KEEP_WORK" == true ]]; then
      log "Work directory kept at: $WORK_DIR"
    else
      rm -rf -- "$WORK_DIR"
    fi
  fi
}

prepare() {
  [[ -n "$INPUT_DIR" ]] || { usage >&2; exit 1; }
  [[ -d "$INPUT_DIR" ]] || die "input directory not found: $INPUT_DIR"
  [[ "$LAYERS" =~ ^[1-9][0-9]*$ ]] || die "--layers must be a positive integer"
  [[ -x "$WDROP_BIN" ]] || die "wdrop executable not found or not executable: $WDROP_BIN (run make in src/)"

  WDROP_BIN="$(cd "$(dirname "$WDROP_BIN")" && pwd)/$(basename "$WDROP_BIN")"
  INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"

  if [[ "$PREPROCESS" == true ]]; then
    [[ -f "$PREPROCESSOR" ]] || die "preprocessor not found: $PREPROCESSOR"
    if [[ -f "$PROJECT_ROOT/.venv/bin/activate" ]]; then
      # shellcheck disable=SC1091
      source "$PROJECT_ROOT/.venv/bin/activate"
    fi
  fi

  if [[ -z "$OUTPUT_CSV" ]]; then
    mkdir -p "$SCRIPT_DIR/results"
    OUTPUT_CSV="$SCRIPT_DIR/results/$(date -u +%Y%m%dT%H%M%SZ).csv"
  else
    mkdir -p "$(dirname "$OUTPUT_CSV")"
  fi

  if [[ -z "$WORK_DIR" ]]; then
    WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/wdrop-layer-runtime.XXXXXX")"
  else
    mkdir -p "$WORK_DIR"
    WORK_DIR="$(cd "$WORK_DIR" && pwd)"
  fi
}

# Runs wdrop once on in.pdb inside RUN_DIR (short relative names: wdrop limits
# file names to 100 chars). Prints the runtime in ms; on success in.pdb is replaced
# by wdrop's output so the next call continues from it.
run_layer() {
  local RUN_DIR="$1"
  local LAYER="$2"
  local OUT="in_1WAT.pdb"
  local LOG_FILE RUNTIME_MS

  (
    cd "$RUN_DIR"
    "$WDROP_BIN" --file in.pdb --sigma "$SIGMA_P" --weed-dist "$WEED_DIST" --layers 1 \
      > "layer${LAYER}.stdout" 2> "layer${LAYER}.stderr"
  ) || return 1

  [[ -f "$RUN_DIR/$OUT" ]] || return 1

  # wdrop writes its log as "<output>-wdrop-YYYY-MM-DD HH:MM:SS"
  LOG_FILE="$(find "$RUN_DIR" -maxdepth 1 -name "${OUT}-wdrop-*" -print -quit)"
  [[ -n "$LOG_FILE" ]] || return 1
  RUNTIME_MS="$(awk '/^runtime:/ {print $2; exit}' "$LOG_FILE")"
  [[ -n "$RUNTIME_MS" ]] || return 1

  mv -- "$LOG_FILE" "$RUN_DIR/layer${LAYER}.log"
  cp -- "$RUN_DIR/in.pdb" "$RUN_DIR/layer$((LAYER - 1)).pdb"
  mv -- "$RUN_DIR/$OUT" "$RUN_DIR/in.pdb"
  printf '%s' "$RUNTIME_MS"
}

main() {
  parse_args "$@"
  prepare
  trap cleanup EXIT

  local PDBS=()
  mapfile -t PDBS < <(find "$INPUT_DIR" -maxdepth 1 -type f -iname '*.pdb' | sort)
  (( ${#PDBS[@]} > 0 )) || die "no .pdb files in $INPUT_DIR"

  local HEADER="pdb_name"
  for ((L = 1; L <= LAYERS; L++)); do HEADER+=",layer${L}_ms"; done
  printf '%s\n' "$HEADER" > "$OUTPUT_CSV"

  log "Input directory : $INPUT_DIR (${#PDBS[@]} PDB files)"
  log "Layers          : $LAYERS (wdrop --layers 1 each, output chained)"
  log "Preprocess      : $PREPROCESS"
  log "Work directory  : $WORK_DIR"
  log "Output CSV      : $OUTPUT_CSV"

  local FAILURES=0
  local PDB NAME RUN_DIR ROW MS
  for PDB in "${PDBS[@]}"; do
    NAME="$(basename "$PDB")"
    NAME="${NAME%.*}"
    RUN_DIR="$WORK_DIR/$NAME"
    mkdir -p "$RUN_DIR"
    cp -- "$PDB" "$RUN_DIR/in.pdb"
    ROW="$NAME"

    if [[ "$PREPROCESS" == true ]] \
      && ! python3 "$PREPROCESSOR" --target "$RUN_DIR/in.pdb" > "$RUN_DIR/preprocess.log" 2>&1; then
      log "$NAME: preprocessing FAILED (see $RUN_DIR/preprocess.log)"
      ROW+=",FAILED"
      for ((L = 2; L <= LAYERS; L++)); do ROW+=","; done
      printf '%s\n' "$ROW" >> "$OUTPUT_CSV"
      FAILURES=$((FAILURES + 1))
      continue
    fi

    for ((L = 1; L <= LAYERS; L++)); do
      if MS="$(run_layer "$RUN_DIR" "$L")"; then
        ROW+=",$MS"
      else
        log "$NAME: layer $L FAILED (see $RUN_DIR/layer${L}.stdout / .stderr)"
        ROW+=",FAILED"
        for ((R = L + 1; R <= LAYERS; R++)); do ROW+=","; done
        FAILURES=$((FAILURES + 1))
        KEEP_WORK=true
        break
      fi
    done

    printf '%s\n' "$ROW" >> "$OUTPUT_CSV"
    log "$ROW"
  done

  cat "$OUTPUT_CSV"
  log "Done: ${#PDBS[@]} PDB files, $FAILURES failed. CSV: $OUTPUT_CSV"
  (( FAILURES == 0 ))
}

main "$@"
