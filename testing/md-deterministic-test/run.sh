#!/usr/bin/env bash
# Driver for the MD determinism suite: runs the three cases and collects their
# verdicts into one summary.
#
# Cases:
#   gromacs      gmx mdrun                    (CPU)
#   tinker       Tinker `dynamic`             (CPU)
#   tinker-gpu   Tinker9 `dynamic9`           (GPU; SKIPPED when dynamic9 is absent)
#
# A failing case does not stop the others; the driver exits non-zero at the end if
# any case failed. A skipped GPU case is not a failure.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_ROOT="$SCRIPT_DIR/results"

ENGINES="all"
PASSTHROUGH=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [-e ENGINE] [OPTIONS]

Measure how much repeated, identical MD runs differ.

  -e, --engine ENGINE    gromacs | tinker | tinker-gpu | all   (default: all)

All other options are passed through to the per-case scripts:
  -n, --replicates N     number of identical runs to compare (default: 5)
  -i, --input PATH       input PDB (default: testing/test_pdbs/1PSV_cryst.pdb)
      --md-ns NS         MD length in nanoseconds (default: 0.02 = 20000 steps)
      --arms WHICH       pipeline | hardened | all   (default: pipeline)
      --keep             keep the temporary work directories
  -h, --help             show this help message

MD length: 0.02 ns is 20000 steps at the pipeline's 1 fs timestep. Do not go below
~0.002 ns -- TINKER-DETERMINISM.md section 2 measured divergence as still invisible
at 500 steps, so a shorter run reports a false pass.

Examples:
  $(basename "$0")
  $(basename "$0") -e tinker -n 3 --md-ns 0.005 --arms all
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -e|--engine)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 1; }
      ENGINES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) PASSTHROUGH+=("$1"); shift ;;
  esac
done

case "$ENGINES" in
  gromacs|tinker|tinker-gpu|all) ;;
  *) echo "ERROR: --engine must be gromacs, tinker, tinker-gpu or all (got '$ENGINES')" >&2; exit 1 ;;
esac

# Parallel arrays rather than packed strings: the case scripts live under paths the
# user controls, and --gpu must stay a separate argument from the script path.
CASE_NAMES=()
CASE_SCRIPTS=()
CASE_FLAGS=()

add_case() {
  CASE_NAMES+=("$1")
  CASE_SCRIPTS+=("$2")
  CASE_FLAGS+=("${3:-}")
}

[[ "$ENGINES" == "gromacs"    || "$ENGINES" == "all" ]] && add_case "gromacs (CPU)" "$SCRIPT_DIR/GROMACS/run.sh"
[[ "$ENGINES" == "tinker"     || "$ENGINES" == "all" ]] && add_case "tinker (CPU)"  "$SCRIPT_DIR/TINKER/run.sh"
[[ "$ENGINES" == "tinker-gpu" || "$ENGINES" == "all" ]] && add_case "tinker (GPU)"  "$SCRIPT_DIR/TINKER/run.sh" "--gpu"

mkdir -p -- "$RESULTS_ROOT"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
SUMMARY="$RESULTS_ROOT/summary_${STAMP}.md"
HEADLINES="$(mktemp)"
trap 'rm -f -- "$HEADLINES"' EXIT

# Marks a case boundary in the collected headlines; no verdict line can start
# with it, so it is unambiguous when the summary is rebuilt below.
readonly CASE_MARKER="### case: "

export DETERMINISM_SUMMARY_FILE="$HEADLINES"

{
  printf '# MD determinism — summary\n\n'
  printf -- '- Generated: %s\n' "$(date -u +'%Y-%m-%d %H:%M:%S UTC')"
  printf -- '- Cases: %s\n\n' "$ENGINES"
} > "$SUMMARY"

FAILED=()
for I in "${!CASE_NAMES[@]}"; do
  NAME="${CASE_NAMES[$I]}"
  SCRIPT="${CASE_SCRIPTS[$I]}"
  FLAG="${CASE_FLAGS[$I]}"

  printf '\n===== MD determinism: %s =====\n' "$NAME" >&2
  printf '%s\n' "$CASE_MARKER$NAME" >> "$HEADLINES"

  if ! "$SCRIPT" ${FLAG:+"$FLAG"} ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}; then
    printf '%s\n' "$NAME: FAILED (see output above)" >> "$HEADLINES"
    FAILED+=("$NAME")
  fi
done

{
  printf '## Verdicts\n\n'
  while IFS= read -r LINE; do
    if [[ "$LINE" == "$CASE_MARKER"* ]]; then
      printf '\n**%s**\n\n' "${LINE#"$CASE_MARKER"}"
    else
      printf -- '- %s\n' "$LINE"
    fi
  done < "$HEADLINES"
  printf '\nPer-case reports are in `%s/<timestamp>_<case>/report.md`.\n' "results"
} >> "$SUMMARY"

printf '\n' >&2
cat "$SUMMARY" >&2
printf '\nSummary written to: %s\n' "$SUMMARY" >&2

if (( ${#FAILED[@]} > 0 )); then
  printf 'FAILED case(s): %s\n' "${FAILED[*]}" >&2
  exit 1
fi
