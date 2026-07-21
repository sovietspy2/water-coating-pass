#!/usr/bin/env bash
# Determinism of the GROMACS energy-minimisation steps.
#
# The pipeline minimises in two passes (pipeline/gromacs.sh:203-213), so both are
# measured separately -- if only one of them is stochastic, that is worth knowing:
#   em  steepest descent   gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log
#   cg  conjugate gradient gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log
#
# grompp runs once per stage, so every replicate starts from a byte-identical .tpr.
# The cg stage starts from one canonical after_em.gro (produced by a dedicated prep
# run) so that cg variance is not contaminated by em variance.
#
# Arms:
#   pipeline  gmx mdrun exactly as the pipeline calls it (default threading)
#   hardened  gmx mdrun -reprod -ntmpi 1 -ntomp 1
# -reprod is GROMACS's own switch for run-to-run reproducibility; Tinker has no
# equivalent (TINKER-DETERMINISM.md section 2), which is the asymmetry these suites
# are here to quantify.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../../determinism-common/determinism_lib.sh
source "$SUITE_DIR/../determinism-common/determinism_lib.sh"

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run the GROMACS steepest-descent and conjugate-gradient minimisations N times each
on byte-identical inputs and report how far the results differ.

Options:
$(common_usage_options)

Arms: 'pipeline' = gmx mdrun as the pipeline calls it,
      'hardened' = gmx mdrun -reprod -ntmpi 1 -ntomp 1. Default: all.
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

command -v gmx   >/dev/null 2>&1 || die "'gmx' not found in PATH"
command -v wdrop >/dev/null 2>&1 || die "'wdrop' not found in PATH (run 'cd src && make install')"

# Keep gmx from littering #backup# files; each replicate already gets a clean dir.
export GMX_MAXBACKUP=-1

require_input_pdb
setup_work_dir
trap determinism_cleanup EXIT
activate_venv

RESULTS_DIR="$(make_results_dir "$SUITE_DIR/results" "gromacs")"
REPORT="$RESULTS_DIR/report.md"
CSV="$RESULTS_DIR/results.csv"

# ---------------------------------------------------------------------------
# Shared fixture: topology, both .tpr files, and a canonical minimised structure.
# ---------------------------------------------------------------------------

FIXTURE_PDB="$(build_fixture_pdb)"

PREP_DIR="$WORK_DIR/prep"
mkdir -p -- "$PREP_DIR"

log "Prep 1/4: pdb2gmx + editconf (topology and pseudo-PBC box)"
gmx_build_topology "$PREP_DIR" "$FIXTURE_PDB"

log "Prep 2/4: grompp for steepest descent"
gmx_write_st_mdp "$PREP_DIR/gromacs-st.mdp"
gmx_write_cg_mdp "$PREP_DIR/gromacs-cg.mdp"
( cd "$PREP_DIR" && gmx grompp -v -f gromacs-st.mdp -c conf_largebox.gro -r conf_largebox.gro \
    -o em -p topol.top -maxwarn 2 ) >"$PREP_DIR/grompp_em.log" 2>&1 \
  || { tail -40 "$PREP_DIR/grompp_em.log" >&2; die "gmx grompp (em) failed"; }

# One canonical em run, so the cg replicates all start from the same structure.
log "Prep 3/4: canonical steepest-descent run (input for the cg stage)"
( cd "$PREP_DIR" && gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log ) \
  >"$PREP_DIR/mdrun_em.log" 2>&1 \
  || { tail -40 "$PREP_DIR/mdrun_em.log" >&2; die "canonical gmx mdrun (em) failed"; }

log "Prep 4/4: grompp for conjugate gradient"
( cd "$PREP_DIR" && gmx grompp -v -f gromacs-cg.mdp -c after_em.gro -r after_em.gro \
    -o cg -p topol.top ) >"$PREP_DIR/grompp_cg.log" 2>&1 \
  || { tail -40 "$PREP_DIR/grompp_cg.log" >&2; die "gmx grompp (cg) failed"; }

TOTAL_ATOMS="$(awk 'NR==2 {print $1; exit}' "$PREP_DIR/conf_largebox.gro")"
ATOM_LABELS="$PREP_DIR/atom_labels.txt"
gmx_write_atom_labels "$PREP_DIR/conf.gro" "$ATOM_LABELS"
log "System: $TOTAL_ATOMS atoms"

write_report_header "$REPORT" "GROMACS minimize — determinism report" \
  "Stage 1 under test: \`gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log\` (steepest descent)" \
  "Stage 2 under test: \`gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log\` (conjugate gradient)" \
  "System: $TOTAL_ATOMS atoms" \
  "Comparison precision: trajectories converted to GROMOS96 (\`%15.9f\` nm); \`.gro\` output is only \`%8.3f\` nm and is not used"

# ---------------------------------------------------------------------------
# Replicates
# ---------------------------------------------------------------------------

# $1 = stage (em|cg), $2 = arm
run_stage_arm() {
  local STAGE="$1"
  local ARM="$2"
  local ARM_DIR="$WORK_DIR/${STAGE}_${ARM}"
  local OUTPUTS=()
  local LOGS=()
  local TIMINGS=()
  local REP REP_DIR
  local ARM_FLAGS=()
  local TIME_PREFIX=()

  mapfile -t ARM_FLAGS < <(gmx_mdrun_arm_flags "$ARM")

  log "=== stage '$STAGE' arm '$ARM': $REPLICATES replicate(s) ==="

  for ((REP = 1; REP <= REPLICATES; REP++)); do
    REP_DIR="$ARM_DIR/rep_$(printf '%02d' "$REP")"
    mkdir -p -- "$REP_DIR"
    mapfile -t TIME_PREFIX < <(timing_prefix "$REP_DIR/cputime.txt")
    cp -- "$PREP_DIR/$STAGE.tpr" "$REP_DIR/$STAGE.tpr"

    log "stage=$STAGE arm=$ARM replicate $REP/$REPLICATES"
    if ! ( cd "$REP_DIR" && ${TIME_PREFIX[@]+"${TIME_PREFIX[@]}"} \
             gmx mdrun -v -s "$STAGE" -o "$STAGE.trr" \
             -c "after_$STAGE.gro" -g "$STAGE.log" \
             ${ARM_FLAGS[@]+"${ARM_FLAGS[@]}"} ) >"$REP_DIR/mdrun.stdout" 2>&1; then
      RUN_FAILED=true
      tail -30 "$REP_DIR/mdrun.stdout" >&2 || true
      die "gmx mdrun ($STAGE) failed on arm=$ARM replicate $REP (see $REP_DIR)"
    fi

    # .gro resolves to 0.01 A, far too coarse to see early divergence -- convert the
    # full-precision trajectory instead.
    if ! ( cd "$REP_DIR" && gmx trjconv -f "$STAGE.trr" -s "$STAGE.tpr" -o final.g96 <<EOF
0
EOF
    ) >"$REP_DIR/trjconv.log" 2>&1; then
      RUN_FAILED=true
      tail -30 "$REP_DIR/trjconv.log" >&2 || true
      die "gmx trjconv failed on arm=$ARM replicate $REP (see $REP_DIR)"
    fi

    OUTPUTS+=("$REP_DIR/final.g96")
    LOGS+=("$REP_DIR/$STAGE.log")
    # Plain `[[ ... ]] && ...` as the last statement in the loop body would return 1
    # and abort the loop under `set -e` whenever the file is absent.
    if [[ -f "$REP_DIR/cputime.txt" ]]; then
      TIMINGS+=("$REP_DIR/cputime.txt")
    fi
  done

  report_cpu_utilisation "$REPORT" "$STAGE / arm=$ARM" "gmx mdrun" ${TIMINGS[@]+"${TIMINGS[@]}"}

  local LABEL="GROMACS minimize ($STAGE) / arm=$ARM"
  local HEADLINE
  HEADLINE="$("$COMPARE_PY" --format g96 --label "$LABEL" \
      --labels "$ATOM_LABELS" --out-md "$REPORT" --out-csv "$CSV" \
      --out-json "$RESULTS_DIR/${STAGE}_${ARM}_structures.json" "${OUTPUTS[@]}")"
  emit_headline "$HEADLINE"

  # A stage that converges immediately did no work, so "deterministic" is vacuous.
  # This happens to cg whenever em already drove Fmax below the cg emtol (750).
  if grep -qE "converged to \S+ < \S+ in [01] steps" "${LOGS[0]}"; then
    log "NOTE: $STAGE converged in <=1 step — its verdict carries no information for this input"
    printf '> NOTE: this stage converged in <=1 step, so it did essentially no work on this\n' >> "$REPORT"
    printf '> input and its determinism verdict is vacuous. On a small system the steepest\n' >> "$REPORT"
    printf '> descent pass already satisfies the cg `emtol`.\n\n' >> "$REPORT"
  fi

  "$PARSE_SCALARS_PY" --kind gromacs-minimize --label "$LABEL" \
    --out-md "$REPORT" --out-json "$RESULTS_DIR/${STAGE}_${ARM}_scalars.json" "${LOGS[@]}"
}

ARMS="$(resolve_arms all)"
for STAGE in em cg; do
  while read -r ARM; do
    [[ -n "$ARM" ]] || continue
    run_stage_arm "$STAGE" "$ARM"
  done <<< "$ARMS"
done

log "Report: $REPORT"
log "CSV:    $CSV"
