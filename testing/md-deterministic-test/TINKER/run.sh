#!/usr/bin/env bash
# Determinism of the Tinker molecular-dynamics step — CPU and GPU cases.
#
# Under test: the `dynamic` invocation from pipeline/tinker.sh:195 --
#     dynamic md_start.xyz -k md.key
#     stdin: <steps> <dt_fs> <save_interval_ps> 2 300      (NVT at 300 K)
#
# Every replicate starts from ONE canonical minimised structure, produced by a single
# prep-stage `minimize` run with `OPENMP-THREADS 1`. Minimize variance therefore
# cannot leak into the MD numbers; what this suite reports is the MD step alone.
#
# Cases:
#   CPU (default)  `dynamic`
#   GPU (--gpu)    `dynamic9`, with the two extra md.key lines pipeline/tinker.sh:169
#                  injects for Tinker9. Reports SKIPPED (exit 0) when dynamic9 is absent.
#
# Arms:
#   pipeline  md.key exactly as shipped (RANDOMSEED fixed, no thread pin)
#   hardened  the same key + `OPENMP-THREADS 1`
# TINKER-DETERMINISM.md section 4 leaves that one line as the open question for CPU MD;
# this is the harness that answers it. It is a no-op for the GPU case (an OpenMP
# keyword cannot affect forces evaluated on the GPU), so --gpu runs one arm only.
#
# MD length matters: section 2 of that document measured divergence still invisible at
# 500 steps and only Angstrom-scale by 1860, so a short run gives a false pass. The
# default of 0.02 ns = 20000 steps is deliberately well past that point.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../../determinism-common/determinism_lib.sh
source "$SUITE_DIR/../determinism-common/determinism_lib.sh"

readonly DT_FS="1.0" # pipeline/tinker.sh:50
USE_GPU=false

usage() {
  cat <<EOF
Usage: $(basename "$0") [--gpu] [OPTIONS]

Run Tinker \`dynamic\` N times from a byte-identical starting structure and report
how far the trajectories diverge.

Options:
      --gpu              use dynamic9 (Tinker9-GPU) instead of dynamic
$(common_usage_options)

Arms: 'pipeline' = md.key as shipped, 'hardened' = md.key + OPENMP-THREADS 1.
Default: pipeline only (MD is expensive); pass --arms all for both.
--gpu always runs a single arm.
EOF
}

parse_common_args "$@"
set -- "${DETERMINISM_UNPARSED_ARGS[@]+"${DETERMINISM_UNPARSED_ARGS[@]}"}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gpu) USE_GPU=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1" ;;
  esac
done
validate_common_args

awk -v D="$MD_NS" 'BEGIN { exit !(D+0 > 0) }' || die "--md-ns must be greater than 0 for the MD suite (got '$MD_NS')"

if [[ "$USE_GPU" == true ]]; then
  DYNAMIC_CMD="dynamic9"
  CASE_SLUG="tinker_gpu"
  CASE_NAME="TINKER MD (GPU)"
  GPU_FLAG="gpu"
  # Not an error: this box may simply have no GPU build. Report and stand down so
  # `run.sh -e all` still completes and the case is visible in the summary.
  if ! command -v dynamic9 >/dev/null 2>&1; then
    printf 'SKIPPED: %s — dynamic9 not found in PATH (no Tinker9-GPU build on this machine)\n' \
      "$CASE_NAME" >&2
    if [[ -n "${DETERMINISM_SUMMARY_FILE:-}" ]]; then
      printf '%s: SKIPPED — dynamic9 not found in PATH\n' "$CASE_NAME" >> "$DETERMINISM_SUMMARY_FILE"
    fi
    exit 0
  fi
else
  DYNAMIC_CMD="dynamic"
  CASE_SLUG="tinker_cpu"
  CASE_NAME="TINKER MD (CPU)"
  GPU_FLAG=""
  command -v dynamic >/dev/null 2>&1 || die "Tinker 'dynamic' not found in PATH"
fi

command -v minimize >/dev/null 2>&1 || die "Tinker 'minimize' not found in PATH"
command -v pdbxyz   >/dev/null 2>&1 || die "Tinker 'pdbxyz' not found in PATH"
command -v wdrop    >/dev/null 2>&1 || die "'wdrop' not found in PATH (run 'cd src && make install')"

require_input_pdb
setup_work_dir
trap determinism_cleanup EXIT
activate_venv

RESULTS_DIR="$(make_results_dir "$SUITE_DIR/results" "$CASE_SLUG")"
REPORT="$RESULTS_DIR/report.md"
CSV="$RESULTS_DIR/results.csv"

# ---------------------------------------------------------------------------
# Shared fixture: structure -> XYZ -> ONE canonical minimisation.
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

log "Prep: canonical minimize (thread-pinned, so the MD start point is fixed)"
cp -- "$SYS_XYZ" "$PREP_DIR/sys.xyz"
tinker_write_minimize_key "$PREP_DIR/minimize.key" "$PROTEIN_ATOMS" "pipeline"
( cd "$PREP_DIR" && minimize sys.xyz -k minimize.key >minimize.log 2>&1 <<EOF
0.01
EOF
) || { tail -30 "$PREP_DIR/minimize.log" >&2; die "canonical minimize failed"; }
[[ -f "$PREP_DIR/sys.xyz_2" ]] || die "canonical minimize produced no sys.xyz_2"

MD_START="$PREP_DIR/md_start.xyz"
cp -- "$PREP_DIR/sys.xyz_2" "$MD_START"

# Dynamics parameters — pipeline/tinker.sh:174-185
N_STEPS="$(ns_to_steps "$MD_NS" "$DT_FS")"
SAVE_EVERY="$(save_every_steps "$N_STEPS" "$TARGET_FRAMES")"
SAVE_INTERVAL_PS="$(awk -v STEPS="$SAVE_EVERY" -v DT="$DT_FS" 'BEGIN { printf "%.6f\n", (STEPS * DT) / 1000.0 }')"
EXPECTED_FRAMES=$(( (N_STEPS + SAVE_EVERY - 1) / SAVE_EVERY ))

log "MD: $MD_NS ns = $N_STEPS steps at $DT_FS fs; saving every $SAVE_EVERY steps (~$EXPECTED_FRAMES frames)"

write_report_header "$REPORT" "$CASE_NAME — determinism report" \
  "Command under test: \`$DYNAMIC_CMD md_start.xyz -k md.key\` (NVT, 300 K)" \
  "MD length: $MD_NS ns = $N_STEPS steps at $DT_FS fs, ~$EXPECTED_FRAMES saved frames" \
  "System: $TOTAL_ATOMS atoms, $PROTEIN_ATOMS restrained protein atoms" \
  "Start structure: one canonical thread-pinned \`minimize\` result, shared by every replicate" \
  "Comparison precision: \`.arc\` is \`%12.6f\` A; \`.dyn\` (velocities, 16 digits) is hashed separately"

# ---------------------------------------------------------------------------
# Replicates
# ---------------------------------------------------------------------------

run_arm() {
  local ARM="$1"
  local ARM_DIR="$WORK_DIR/md_$ARM"
  local OUTPUTS=()
  local LOGS=()
  local DYN_HASHES=()
  local TIMINGS=()
  local REP REP_DIR
  local TIME_PREFIX=()

  log "=== arm '$ARM': $REPLICATES replicate(s) ==="

  for ((REP = 1; REP <= REPLICATES; REP++)); do
    REP_DIR="$ARM_DIR/rep_$(printf '%02d' "$REP")"
    mkdir -p -- "$REP_DIR"

    cp -- "$MD_START" "$REP_DIR/md_start.xyz"
    cp -- "$PIPELINE_DIR/amber99.prm" "$REP_DIR/amber99.prm"
    tinker_write_md_key "$REP_DIR/md.key" "$PROTEIN_ATOMS" "$ARM" "$GPU_FLAG"

    mapfile -t TIME_PREFIX < <(timing_prefix "$REP_DIR/cputime.txt")

    log "arm=$ARM replicate $REP/$REPLICATES ($N_STEPS steps)"
    if ! ( cd "$REP_DIR" && ${TIME_PREFIX[@]+"${TIME_PREFIX[@]}"} \
             "$DYNAMIC_CMD" md_start.xyz -k md.key >dynamic.log 2>&1 <<EOF
${N_STEPS}
${DT_FS}
${SAVE_INTERVAL_PS}
2
300
EOF
    ); then
      RUN_FAILED=true
      tail -30 "$REP_DIR/dynamic.log" >&2 || true
      die "$DYNAMIC_CMD failed on arm=$ARM replicate $REP (see $REP_DIR)"
    fi

    [[ -f "$REP_DIR/md_start.arc" ]] \
      || { RUN_FAILED=true; die "$DYNAMIC_CMD produced no md_start.arc in $REP_DIR"; }

    OUTPUTS+=("$REP_DIR/md_start.arc")
    LOGS+=("$REP_DIR/dynamic.log")
    # Plain `[[ ... ]] && ...` as the last statement in the loop body would return 1
    # and abort the loop under `set -e` whenever the file is absent.
    if [[ -f "$REP_DIR/cputime.txt" ]]; then
      TIMINGS+=("$REP_DIR/cputime.txt")
    fi
    if [[ -f "$REP_DIR/md_start.dyn" ]]; then
      DYN_HASHES+=("$(sha256sum "$REP_DIR/md_start.dyn" | awk '{print $1}')")
    fi
  done

  report_cpu_utilisation "$REPORT" "arm=$ARM" "$DYNAMIC_CMD" ${TIMINGS[@]+"${TIMINGS[@]}"}

  local LABEL="$CASE_NAME / arm=$ARM"
  local HEADLINE
  HEADLINE="$("$COMPARE_PY" --format tinker-xyz --trajectory --label "$LABEL" \
      --out-md "$REPORT" --out-csv "$CSV" --out-json "$RESULTS_DIR/${ARM}_structures.json" \
      "${OUTPUTS[@]}")"
  emit_headline "$HEADLINE"

  # The .arc archive keeps only 6 decimals. The restart file keeps 16, so it detects
  # divergence that is still invisible in the trajectory (TINKER-DETERMINISM.md s.2).
  if (( ${#DYN_HASHES[@]} > 0 )); then
    local DISTINCT
    DISTINCT="$(printf '%s\n' "${DYN_HASHES[@]}" | sort -u | wc -l)"
    {
      printf '**Restart file (`md_start.dyn`, 16 significant digits) — %s distinct of %s runs.**\n\n' \
        "$DISTINCT" "${#DYN_HASHES[@]}"
      printf 'This is the most sensitive detector available: velocities and accelerations are\n'
      printf 'written at full precision, whereas the `.arc` trajectory keeps only 6 decimals.\n\n'
    } >> "$REPORT"
    emit_headline "$LABEL: .dyn restart files — $DISTINCT/${#DYN_HASHES[@]} distinct"
  fi

  "$PARSE_SCALARS_PY" --kind tinker-dynamic --label "$LABEL" \
    --out-md "$REPORT" --out-json "$RESULTS_DIR/${ARM}_scalars.json" "${LOGS[@]}"
}

if [[ "$USE_GPU" == true ]]; then
  if [[ -n "$ARMS_SELECTOR" && "$ARMS_SELECTOR" != "pipeline" ]]; then
    log "NOTE: ignoring --arms $ARMS_SELECTOR — OPENMP-THREADS is an OpenMP/CPU keyword and"
    log "      cannot affect forces evaluated on the GPU, so --gpu has a single arm."
  fi
  ARMS="$(resolve_arms pipeline single)"
else
  ARMS="$(resolve_arms pipeline)"
fi

while read -r ARM; do
  [[ -n "$ARM" ]] || continue
  run_arm "$ARM"
done <<< "$ARMS"

log "Report: $REPORT"
log "CSV:    $CSV"
