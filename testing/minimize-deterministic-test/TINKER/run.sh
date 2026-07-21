#!/usr/bin/env bash
# Determinism of the Tinker energy-minimisation step.
#
# Under test: the `minimize` invocation from pipeline/tinker.sh:140 --
#     minimize sys.xyz -k minimize.key   (gradient criterion 0.01 on stdin)
# The binary is called directly rather than through tinker.sh so that every input
# byte -- .xyz, .key, amber99.prm -- is identical across replicates and the only
# thing that can vary is the engine itself.
#
# Arms:
#   pipeline  minimize.key WITH `OPENMP-THREADS 1` (what the pipeline ships today)
#   hardened  the same key WITHOUT that line, i.e. default thread count
# Note the direction: for minimize the pipeline is already the deterministic one,
# so the "hardened" arm here is the control that is expected to FAIL. It exists to
# prove the line is still doing work on this machine -- see TINKER-DETERMINISM.md
# section 4 and section 5.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../../determinism-common/determinism_lib.sh
source "$SUITE_DIR/../determinism-common/determinism_lib.sh"

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run Tinker \`minimize\` N times on byte-identical inputs and report how far the
results differ.

Options:
$(common_usage_options)

Arms: 'pipeline' = minimize.key with OPENMP-THREADS 1 (as shipped),
      'hardened' = without it (default threads). Default: all.
EOF
}

parse_common_args "$@"
set -- "${DETERMINISM_UNPARSED_ARGS[@]+"${DETERMINISM_UNPARSED_ARGS[@]}"}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1" ;;
  esac
done
validate_common_args

command -v minimize >/dev/null 2>&1 || die "Tinker 'minimize' not found in PATH"
command -v pdbxyz   >/dev/null 2>&1 || die "Tinker 'pdbxyz' not found in PATH"
command -v wdrop    >/dev/null 2>&1 || die "'wdrop' not found in PATH (run 'cd src && make install')"

require_input_pdb
setup_work_dir
trap determinism_cleanup EXIT
activate_venv

RESULTS_DIR="$(make_results_dir "$SUITE_DIR/results" "tinker")"
REPORT="$RESULTS_DIR/report.md"
CSV="$RESULTS_DIR/results.csv"

# ---------------------------------------------------------------------------
# Shared fixture: one prepared structure, converted to XYZ exactly once.
# ---------------------------------------------------------------------------

FIXTURE_PDB="$(build_fixture_pdb)"

PREP_DIR="$WORK_DIR/prep"
mkdir -p -- "$PREP_DIR"
SYS_XYZ="$(tinker_build_xyz "$PREP_DIR" "$FIXTURE_PDB")"
PROTEIN_ATOMS="$(tinker_protein_atom_count "$SYS_XYZ")"
[[ -n "$PROTEIN_ATOMS" && "$PROTEIN_ATOMS" -ge 1 ]] \
  || die "could not determine the protein atom count from $SYS_XYZ (no OW/HW atoms found)"

TOTAL_ATOMS="$(awk 'NR==1 {print $1; exit}' "$SYS_XYZ")"
log "System: $TOTAL_ATOMS atoms, $PROTEIN_ATOMS restrained protein atoms"

write_report_header "$REPORT" "Tinker minimize — determinism report" \
  "Command under test: \`minimize sys.xyz -k minimize.key\` (gradient criterion 0.01)" \
  "System: $TOTAL_ATOMS atoms, $PROTEIN_ATOMS restrained protein atoms" \
  "Comparison precision: Tinker \`.xyz\` is \`%12.6f\` A"

# ---------------------------------------------------------------------------
# Replicates
# ---------------------------------------------------------------------------

run_arm() {
  local ARM="$1"
  local ARM_DIR="$WORK_DIR/minimize_$ARM"
  local OUTPUTS=()
  local LOGS=()
  local TIMINGS=()
  local REP REP_DIR
  local TIME_PREFIX=()

  log "=== arm '$ARM': $REPLICATES replicate(s) ==="

  for ((REP = 1; REP <= REPLICATES; REP++)); do
    REP_DIR="$ARM_DIR/rep_$(printf '%02d' "$REP")"
    mkdir -p -- "$REP_DIR"

    cp -- "$SYS_XYZ" "$REP_DIR/sys.xyz"
    cp -- "$PIPELINE_DIR/amber99.prm" "$REP_DIR/amber99.prm"
    tinker_write_minimize_key "$REP_DIR/minimize.key" "$PROTEIN_ATOMS" "$ARM"

    mapfile -t TIME_PREFIX < <(timing_prefix "$REP_DIR/cputime.txt")

    log "arm=$ARM replicate $REP/$REPLICATES"
    if ! ( cd "$REP_DIR" && ${TIME_PREFIX[@]+"${TIME_PREFIX[@]}"} \
             minimize sys.xyz -k minimize.key >minimize.log 2>&1 <<EOF
0.01
EOF
    ); then
      RUN_FAILED=true
      tail -30 "$REP_DIR/minimize.log" >&2 || true
      die "minimize failed on arm=$ARM replicate $REP (see $REP_DIR)"
    fi

    [[ -f "$REP_DIR/sys.xyz_2" ]] || { RUN_FAILED=true; die "minimize produced no sys.xyz_2 in $REP_DIR"; }

    OUTPUTS+=("$REP_DIR/sys.xyz_2")
    LOGS+=("$REP_DIR/minimize.log")
    # Plain `[[ ... ]] && ...` as the last statement in the loop body would return 1
    # and abort the loop under `set -e` whenever the file is absent.
    if [[ -f "$REP_DIR/cputime.txt" ]]; then
      TIMINGS+=("$REP_DIR/cputime.txt")
    fi
  done

  report_cpu_utilisation "$REPORT" "arm=$ARM" "minimize" ${TIMINGS[@]+"${TIMINGS[@]}"}

  local LABEL="TINKER minimize / arm=$ARM"
  local HEADLINE
  HEADLINE="$("$COMPARE_PY" --format tinker-xyz --label "$LABEL" \
      --out-md "$REPORT" --out-csv "$CSV" --out-json "$RESULTS_DIR/${ARM}_structures.json" \
      "${OUTPUTS[@]}")"
  emit_headline "$HEADLINE"

  "$PARSE_SCALARS_PY" --kind tinker-minimize --label "$LABEL" \
    --out-md "$REPORT" --out-json "$RESULTS_DIR/${ARM}_scalars.json" "${LOGS[@]}"
}

ARMS="$(resolve_arms all)"
while read -r ARM; do
  [[ -n "$ARM" ]] || continue
  run_arm "$ARM"
done <<< "$ARMS"

log "Report: $REPORT"
log "CSV:    $CSV"
