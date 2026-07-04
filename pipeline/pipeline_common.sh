#!/bin/bash
set -euo pipefail

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

validate_script_dir_not_input_dir() {
  local input_path="$1"
  local input_dir_candidate
  local input_dir_real
  local SCRIPT_DIR_VAR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

  input_dir_candidate="$(dirname -- "$input_path")"

  if [[ -d "$input_dir_candidate" ]]; then
    input_dir_real="$(cd -P "$input_dir_candidate" && pwd)"

    if [[ "$SCRIPT_DIR_VAR" == "$input_dir_real" ]]; then
      log "Error: wdrop.sh cannot be located in the same directory as INPUT_PDB." >&2
      log "Move the input PDB to a different folder!" >&2
      exit 1
    fi
  fi
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

# MD (and the MobyWat step it feeds) runs only when a positive duration is set.
# Shared by the MM backends (gromacs.sh, tinker.sh) to gate the MD/MobyWat steps.
md_enabled() { awk -v D="${1:-0}" 'BEGIN { exit !(D+0 > 0) }'; }

run_mm_step() {
  local mode="$1"
  local input_pdb="$2"
  local script_dir="$3"
  local md_duration="${4:-1}" # in nanosec, default 1; >0 runs MD + MobyWat, 0 skips both
  local REFERENCE_PDB_COM="${5}"

  case "$mode" in
    gromacs)
      "$script_dir/gromacs.sh" "$input_pdb" "$md_duration" 1000 "$REFERENCE_PDB_COM"
      ;;
    tinker)
      "$script_dir/tinker.sh" "$input_pdb" "$md_duration" 1000 "$REFERENCE_PDB_COM"
      ;;
    *)
      echo "Error: unsupported mode: $mode (use gromacs or tinker)" >&2
      return 1
      ;;
  esac
}