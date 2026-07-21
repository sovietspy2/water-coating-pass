#!/usr/bin/env bash
# Determinism of the GROMACS molecular-dynamics step.
#
# Under test: the mdrun call from pipeline/gromacs.sh:244 --
#     gmx mdrun -v -s md -o md.trr -c after_md.gro -g md.log
#
# `md.tpr` is built ONCE during prep, so the `gen_vel` starting velocities
# (gen_seed 28480426) are byte-identical in every replicate -- the tpr is the whole
# input state. Prep also runs the steepest-descent and conjugate-gradient passes once
# each, so every replicate starts MD from the same minimised structure and minimize
# variance cannot leak into the MD numbers.
#
# Arms:
#   pipeline  gmx mdrun exactly as the pipeline calls it (default threading)
#   hardened  gmx mdrun -reprod -ntmpi 1 -ntomp 1
#
# MD length: the default 0.02 ns = 20000 steps at the pipeline's 1 fs timestep,
# matching the Tinker case so the two engines are compared over the same simulated
# time. TINKER-DETERMINISM.md section 2 shows shorter runs give false passes.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../../determinism-common/determinism_lib.sh
source "$SUITE_DIR/../determinism-common/determinism_lib.sh"

readonly DT_PS="0.001" # pipeline/gromacs.sh:51

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run GROMACS MD N times from a byte-identical .tpr and report how far the
trajectories diverge.

Options:
$(common_usage_options)

Arms: 'pipeline' = gmx mdrun as the pipeline calls it,
      'hardened' = gmx mdrun -reprod -ntmpi 1 -ntomp 1.
Default: pipeline only (MD is expensive); pass --arms all for both.
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

awk -v D="$MD_NS" 'BEGIN { exit !(D+0 > 0) }' || die "--md-ns must be greater than 0 for the MD suite (got '$MD_NS')"

command -v gmx   >/dev/null 2>&1 || die "'gmx' not found in PATH"
command -v wdrop >/dev/null 2>&1 || die "'wdrop' not found in PATH (run 'cd src && make install')"

export GMX_MAXBACKUP=-1

require_input_pdb
setup_work_dir
trap determinism_cleanup EXIT
activate_venv

RESULTS_DIR="$(make_results_dir "$SUITE_DIR/results" "gromacs")"
REPORT="$RESULTS_DIR/report.md"
CSV="$RESULTS_DIR/results.csv"

# ---------------------------------------------------------------------------
# Shared fixture: topology, one minimisation pass, one md.tpr.
# ---------------------------------------------------------------------------

FIXTURE_PDB="$(build_fixture_pdb)"

PREP_DIR="$WORK_DIR/prep"
mkdir -p -- "$PREP_DIR"

log "Prep 1/4: pdb2gmx + editconf (topology and pseudo-PBC box)"
gmx_build_topology "$PREP_DIR" "$FIXTURE_PDB"

log "Prep 2/4: steepest descent"
gmx_write_st_mdp "$PREP_DIR/gromacs-st.mdp"
( cd "$PREP_DIR" \
    && gmx grompp -v -f gromacs-st.mdp -c conf_largebox.gro -r conf_largebox.gro \
         -o em -p topol.top -maxwarn 2 \
    && gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log ) \
  >"$PREP_DIR/em.stdout" 2>&1 \
  || { tail -40 "$PREP_DIR/em.stdout" >&2; die "steepest-descent prep failed"; }

log "Prep 3/4: conjugate gradient"
gmx_write_cg_mdp "$PREP_DIR/gromacs-cg.mdp"
( cd "$PREP_DIR" \
    && gmx grompp -v -f gromacs-cg.mdp -c after_em.gro -r after_em.gro -o cg -p topol.top \
    && gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log ) \
  >"$PREP_DIR/cg.stdout" 2>&1 \
  || { tail -40 "$PREP_DIR/cg.stdout" >&2; die "conjugate-gradient prep failed"; }

# MD parameters — pipeline/gromacs.sh:220-239
N_STEPS="$(ns_to_steps "$MD_NS" "$(awk -v P="$DT_PS" 'BEGIN{printf "%.6f", P*1000.0}')")"
SAVE_EVERY="$(save_every_steps "$N_STEPS" "$TARGET_FRAMES")"
EXPECTED_FRAMES=$(( (N_STEPS + SAVE_EVERY - 1) / SAVE_EVERY ))

log "Prep 4/4: grompp for MD ($MD_NS ns = $N_STEPS steps, saving every $SAVE_EVERY)"
gmx_write_md_mdp "$PREP_DIR/gromacs-md.mdp" "$N_STEPS" "$SAVE_EVERY" "$DT_PS"
( cd "$PREP_DIR" && gmx grompp -f gromacs-md.mdp -o md -c after_cg.gro -r after_cg.gro \
    -p topol.top -maxwarn 1 ) >"$PREP_DIR/grompp_md.log" 2>&1 \
  || { tail -40 "$PREP_DIR/grompp_md.log" >&2; die "gmx grompp (md) failed"; }

TOTAL_ATOMS="$(awk 'NR==2 {print $1; exit}' "$PREP_DIR/conf_largebox.gro")"
ATOM_LABELS="$PREP_DIR/atom_labels.txt"
gmx_write_atom_labels "$PREP_DIR/conf.gro" "$ATOM_LABELS"
log "System: $TOTAL_ATOMS atoms"

write_report_header "$REPORT" "GROMACS MD — determinism report" \
  "Command under test: \`gmx mdrun -v -s md -o md.trr -c after_md.gro -g md.log\`" \
  "MD length: $MD_NS ns = $N_STEPS steps at $DT_PS ps, ~$EXPECTED_FRAMES saved frames" \
  "System: $TOTAL_ATOMS atoms" \
  "Start state: one canonical \`md.tpr\` (gen_seed 28480426), so the initial velocities are identical in every replicate" \
  "Comparison precision: trajectories converted to GROMOS96 (\`%15.9f\` nm)"

# ---------------------------------------------------------------------------
# Replicates
# ---------------------------------------------------------------------------

run_arm() {
  local ARM="$1"
  local ARM_DIR="$WORK_DIR/md_$ARM"
  local OUTPUTS=()
  local LOGS=()
  local TIMINGS=()
  local REP REP_DIR
  local ARM_FLAGS=()
  local TIME_PREFIX=()

  mapfile -t ARM_FLAGS < <(gmx_mdrun_arm_flags "$ARM")

  log "=== arm '$ARM': $REPLICATES replicate(s) ==="

  for ((REP = 1; REP <= REPLICATES; REP++)); do
    REP_DIR="$ARM_DIR/rep_$(printf '%02d' "$REP")"
    mkdir -p -- "$REP_DIR"
    mapfile -t TIME_PREFIX < <(timing_prefix "$REP_DIR/cputime.txt")
    cp -- "$PREP_DIR/md.tpr" "$REP_DIR/md.tpr"

    log "arm=$ARM replicate $REP/$REPLICATES ($N_STEPS steps)"
    if ! ( cd "$REP_DIR" && ${TIME_PREFIX[@]+"${TIME_PREFIX[@]}"} \
             gmx mdrun -v -s md -o md.trr -c after_md.gro -g md.log \
             ${ARM_FLAGS[@]+"${ARM_FLAGS[@]}"} ) >"$REP_DIR/mdrun.stdout" 2>&1; then
      RUN_FAILED=true
      tail -30 "$REP_DIR/mdrun.stdout" >&2 || true
      die "gmx mdrun (md) failed on arm=$ARM replicate $REP (see $REP_DIR)"
    fi

    # Full trajectory at full precision: the divergence-vs-frame table needs early
    # frames, where the difference is far below what .gro's 0.01 A can represent.
    if ! ( cd "$REP_DIR" && gmx trjconv -f md.trr -s md.tpr -o traj.g96 <<EOF
0
EOF
    ) >"$REP_DIR/trjconv.log" 2>&1; then
      RUN_FAILED=true
      tail -30 "$REP_DIR/trjconv.log" >&2 || true
      die "gmx trjconv failed on arm=$ARM replicate $REP (see $REP_DIR)"
    fi

    OUTPUTS+=("$REP_DIR/traj.g96")
    LOGS+=("$REP_DIR/md.log")
    # Plain `[[ ... ]] && ...` as the last statement in the loop body would return 1
    # and abort the loop under `set -e` whenever the file is absent.
    if [[ -f "$REP_DIR/cputime.txt" ]]; then
      TIMINGS+=("$REP_DIR/cputime.txt")
    fi
  done

  report_cpu_utilisation "$REPORT" "arm=$ARM" "gmx mdrun" ${TIMINGS[@]+"${TIMINGS[@]}"}

  local LABEL="GROMACS MD / arm=$ARM"
  local HEADLINE
  HEADLINE="$("$COMPARE_PY" --format g96 --trajectory --label "$LABEL" \
      --labels "$ATOM_LABELS" --out-md "$REPORT" --out-csv "$CSV" --out-json "$RESULTS_DIR/${ARM}_structures.json" \
      "${OUTPUTS[@]}")"
  emit_headline "$HEADLINE"

  "$PARSE_SCALARS_PY" --kind gromacs-md --label "$LABEL" \
    --out-md "$REPORT" --out-json "$RESULTS_DIR/${ARM}_scalars.json" "${LOGS[@]}"
}

ARMS="$(resolve_arms pipeline)"
while read -r ARM; do
  [[ -n "$ARM" ]] || continue
  run_arm "$ARM"
done <<< "$ARMS"

log "Report: $REPORT"
log "CSV:    $CSV"
