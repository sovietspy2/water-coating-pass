#!/usr/bin/env bash
# Driver for the minimize determinism suite: runs the per-engine cases and
# collects their verdicts into one summary.
#
# Cases:
#   tinker    Tinker `minimize`
#   gromacs   `gmx mdrun` steepest descent + conjugate gradient
#
# A failing case does not stop the others; the driver exits non-zero at the end if
# any case failed, so one broken engine still leaves you the other engine's numbers.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_ROOT="$SCRIPT_DIR/results"

ENGINES="all"
PASSTHROUGH=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [-e ENGINE] [OPTIONS]

Measure how much repeated, identical energy-minimisation runs differ.

  -e, --engine ENGINE    gromacs | tinker | all   (default: all)

All other options are passed through to the per-engine scripts:
  -n, --replicates N     number of identical runs to compare (default: 5)
  -i, --input PATH       input PDB (default: testing/test_pdbs/1PSV_cryst.pdb)
      --arms WHICH       pipeline | hardened | all   (default: all)
      --keep             keep the temporary work directories
  -h, --help             show this help message

Examples:
  $(basename "$0")
  $(basename "$0") -e tinker --arms all -n 3
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
  gromacs|tinker|all) ;;
  *) echo "ERROR: --engine must be gromacs, tinker or all (got '$ENGINES')" >&2; exit 1 ;;
esac

CASES=()
[[ "$ENGINES" == "tinker"  || "$ENGINES" == "all" ]] && CASES+=("tinker:$SCRIPT_DIR/TINKER/run.sh")
[[ "$ENGINES" == "gromacs" || "$ENGINES" == "all" ]] && CASES+=("gromacs:$SCRIPT_DIR/GROMACS/run.sh")

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
  printf '# Minimize determinism — summary\n\n'
  printf -- '- Generated: %s\n' "$(date -u +'%Y-%m-%d %H:%M:%S UTC')"
  printf -- '- Cases: %s\n\n' "$ENGINES"
} > "$SUMMARY"

FAILED=()
for ENTRY in "${CASES[@]}"; do
  NAME="${ENTRY%%:*}"
  SCRIPT="${ENTRY#*:}"

  printf '\n===== minimize determinism: %s =====\n' "$NAME" >&2
  printf '%s\n' "$CASE_MARKER$NAME" >> "$HEADLINES"

  if ! "$SCRIPT" ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}; then
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
