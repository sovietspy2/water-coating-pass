#!/bin/bash
set -euo pipefail

pipeline_script_dir() {
  local src="${BASH_SOURCE[0]}"
  while [[ -h "$src" ]]; do
    local dir
    dir="$(cd -P "$(dirname "$src")" && pwd)"
    src="$(readlink "$src")"
    [[ "$src" != /* ]] && src="$dir/$src"
  done
  cd -P "$(dirname "$src")" && pwd
}

normalize_input_pdb() {
  local input="$1"

  if [[ ! -f "$input" ]]; then
    echo "Error: file not found: $input" >&2
    return 1
  fi

  if [[ "$input" != *.pdb ]]; then
    echo "Error: input must end with .pdb (got: $input)" >&2
    return 1
  fi

  printf '%s\n' "$(cd -P "$(dirname "$input")" && pwd)/$(basename "$input")"
}

make_output_dir() {
  local input_abs="$1"
  local suffix="$2"
  local base_name
  base_name="$(basename "${input_abs%.pdb}")"
  printf '%s/%s%s\n' "$(dirname "$input_abs")" "$base_name" "$suffix"
}

# ============================================================================
# LOGGING FUNCTIONS
# ============================================================================
# These functions manage unified logging across all pipeline scripts.
# Global variable: LOGFILE (set by calling scripts)

setup_logging() {
  local logfile="$1"
  export LOGFILE="$logfile"

  # Initialize log file only if it doesn't exist (don't truncate existing files)
  if [[ ! -f "$LOGFILE" ]]; then
    : > "$LOGFILE"
  fi
}

log() {
  local line
  if [[ -z "${LOGFILE:-}" ]]; then
    while IFS= read -r line; do
      printf '[%(%F %T)T] %s\n' -1 "$line" >&2
    done <<< "$*"
  else
    while IFS= read -r line; do
      printf '[%(%F %T)T] %s\n' -1 "$line"
    done <<< "$*" | tee -a "$LOGFILE"
  fi
}

run_step() {
  log "RUN: $*"
  "$@"
  local rc=$?
  log "DONE (rc=$rc): $*"
  return "$rc"
}

# ============================================================================

run_mm_step() {
  local mode="$1"
  local input_pdb="$2"
  local script_dir="$3"

  case "$mode" in
    gromacs)
      "$script_dir/gromacs.sh" "$input_pdb"
      ;;
    tinker)
      "$script_dir/tinker.sh" "$input_pdb"
      ;;
    *)
      echo "Error: unsupported mode: $mode (use gromacs or tinker)" >&2
      return 1
      ;;
  esac
}