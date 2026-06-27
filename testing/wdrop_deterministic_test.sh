#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

WDROP_BIN="$PROJECT_ROOT/src/wdrop"
INPUT_PDB="$PROJECT_ROOT/src/1PSV.pdb"
N=1000
SIGMA_P="1.8"
WEED_DIST="3.5"
LAYERS="5"
WORK_DIR=""
CLEAN_ON_SUCCESS=false
PROGRESS_EVERY=100

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run wdrop repeatedly with the same input, hash each generated PDB, and fail if
any hash differs from the first run.

Default command under test:
  wdrop --file input.pdb --sigma 1.8 --weed-dist 3.5 --layers 5

Options:
  -n, --runs N           Number of runs (default: 1000)
  -p, --wdrop-bin PATH   wdrop executable path (default: ../src/wdrop)
  -i, --input PATH       input PDB path (default: ../src/1PSV.pdb)
  -w, --work-dir DIR     directory for per-run files (default: mktemp in /tmp)
  --clean                remove work directory after a successful run
  --sigma VALUE          probe radius, default: 1.8
  --weed-dist VALUE      minimal O-O distance, default: 3.5
  --layers N             number of water layers, default: 5
  --progress-every N     print progress every N runs, 0 disables (default: 100)
  -h, --help             show this help message

Examples:
  $(basename "$0")
  $(basename "$0") -n 100
  $(basename "$0") --runs 500 --work-dir ./determinism_runs
EOF
}

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

is_positive_integer() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -n|--runs)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        N="$2"
        shift 2
        ;;
      -p|--wdrop-bin)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        WDROP_BIN="$2"
        shift 2
        ;;
      -i|--input)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        INPUT_PDB="$2"
        shift 2
        ;;
      -w|--work-dir)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        WORK_DIR="$2"
        shift 2
        ;;
      --clean)
        CLEAN_ON_SUCCESS=true
        shift
        ;;
      --sigma)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        SIGMA_P="$2"
        shift 2
        ;;
      --weed-dist)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        WEED_DIST="$2"
        shift 2
        ;;
      --layers)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        LAYERS="$2"
        shift 2
        ;;
      --progress-every)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        PROGRESS_EVERY="$2"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "unknown option: $1"
        ;;
    esac
  done
}

cleanup() {
  local status=$?

  if [[ "$status" -eq 0 && "$CLEAN_ON_SUCCESS" == true && -n "$WORK_DIR" && -d "$WORK_DIR" ]]; then
    rm -rf "$WORK_DIR"
  elif [[ -n "${WORK_DIR:-}" && -d "$WORK_DIR" ]]; then
    log "Work directory kept at: $WORK_DIR"
  fi
}

prepare() {
  is_positive_integer "$N" || die "--runs must be a positive integer"
  is_positive_integer "$LAYERS" || die "--layers must be a positive integer"
  [[ "$PROGRESS_EVERY" =~ ^[0-9]+$ ]] || die "--progress-every must be a non-negative integer"

  [[ -x "$WDROP_BIN" ]] || die "wdrop executable not found or not executable: $WDROP_BIN"
  [[ -f "$INPUT_PDB" ]] || die "input PDB not found: $INPUT_PDB"
  command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required"

  WDROP_BIN="$(cd "$(dirname "$WDROP_BIN")" && pwd)/$(basename "$WDROP_BIN")"
  INPUT_PDB="$(cd "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"

  if [[ -z "$WORK_DIR" ]]; then
    WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/wdrop-determinism.XXXXXX")"
  else
    mkdir -p "$WORK_DIR"
    WORK_DIR="$(cd "$WORK_DIR" && pwd)"
  fi
}

output_name_for_input() {
  local input_name="$1"
  local base="${input_name%.[Pp][Dd][Bb]}"

  if [[ "$LAYERS" == "1" ]]; then
    printf '%s_1WAT.pdb' "$base"
  else
    printf '%s_%sWAT.pdb' "$base" "$LAYERS"
  fi
}

run_once() {
  local run_no="$1"
  local run_dir="$WORK_DIR/run_$(printf '%06d' "$run_no")"
  local input_name="input.pdb"
  local output_name
  output_name="$(output_name_for_input "$input_name")"

  mkdir -p "$run_dir"
  cp "$INPUT_PDB" "$run_dir/$input_name"

  if (
    cd "$run_dir"
    "$WDROP_BIN" --file "$input_name" --sigma "$SIGMA_P" --weed-dist "$WEED_DIST" --layers "$LAYERS" \
      >stdout.txt 2>stderr.txt
  ); then
    :
  else
    local status=$?
    printf 'wdrop failed in run %06d with exit code %d\n' "$run_no" "$status" >&2
    printf 'Run directory: %s\n' "$run_dir" >&2
    if [[ -s "$run_dir/stderr.txt" ]]; then
      printf 'Last stderr lines:\n' >&2
      tail -20 "$run_dir/stderr.txt" >&2
    fi
    return "$status"
  fi

  [[ -f "$run_dir/$output_name" ]] || die "run $run_no did not create expected output: $run_dir/$output_name"

  sha256sum "$run_dir/$output_name" | awk '{print $1}'
}

main() {
  parse_args "$@"
  prepare
  trap cleanup EXIT

  local hashes_file="$WORK_DIR/hashes.tsv"
  local baseline_hash=""
  local mismatches=0

  : > "$hashes_file"

  log "Running $N deterministic checks"
  log "Command: $WDROP_BIN --file input.pdb --sigma $SIGMA_P --weed-dist $WEED_DIST --layers $LAYERS"
  log "Input: $INPUT_PDB"
  log "Work directory: $WORK_DIR"

  for ((run_no = 1; run_no <= N; run_no++)); do
    local hash
    hash="$(run_once "$run_no")"
    printf '%06d\t%s\n' "$run_no" "$hash" >> "$hashes_file"

    if [[ "$run_no" -eq 1 ]]; then
      baseline_hash="$hash"
      log "Baseline hash: $baseline_hash"
    elif [[ "$hash" != "$baseline_hash" ]]; then
      ((mismatches += 1))
      printf 'MISMATCH run=%06d baseline=%s actual=%s dir=%s\n' \
        "$run_no" "$baseline_hash" "$hash" "$WORK_DIR/run_$(printf '%06d' "$run_no")" >&2
    fi

    if (( PROGRESS_EVERY > 0 && run_no % PROGRESS_EVERY == 0 )); then
      log "Progress: $run_no/$N runs completed"
    fi
  done

  if (( mismatches > 0 )); then
    log "FAILED: $mismatches of $N runs differed from the baseline"
    log "Hashes: $hashes_file"
    exit 1
  fi

  log "PASS: all $N runs produced identical generated PDB hashes"
  log "Hashes: $hashes_file"
}

main "$@"
