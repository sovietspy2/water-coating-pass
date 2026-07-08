#!/bin/bash
# ============================================================================
# research.sh — batch benchmark + report generator for the WDROP pipeline.
#
# Processes every .pdb file found in an input directory (the files are provided
# and reviewed manually — this script does NOT download anything). For each PDB
# it runs all four pipeline combinations:
#   tinker default | tinker per_layer | gromacs default | gromacs per_layer
#   (refinement = whether MD runs only on the final layer or after every layer)
# Each combo runs in its OWN working directory (the _l{L}_{refinement} output tags collide
# between engines, and each combo needs its own application.LOG for runtime
# extraction). To run Tinker in GPU mode, export TINKER_GPU=1 before invoking
# this script.
#
# Results are distilled into a single shared Markdown report. For every PDB the
# report shows 4 data areas, each containing:
#   - total runtime (last "Total runtime was N seconds" line of application.LOG)
#   - the final SR value from every O_system_mt*.lst file
#   - the average and maximum C-alpha RMSD from O_system_rmst.txt
#
# Usage:
#   ./pipeline/research.sh <INPUT_DIR> <OUTPUT_DIR>
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

usage() {
  cat >&2 <<EOF
Usage: $(basename "$0") <INPUT_DIR> <OUTPUT_DIR>

  <INPUT_DIR>    directory containing the .pdb files to process (read-only)
  <OUTPUT_DIR>   directory where working files and research_report.md are written
EOF
  exit 1
}

# ----------------------------------------------------------------------------
# Args
# ----------------------------------------------------------------------------
[[ $# -eq 2 ]] || usage
INPUT_DIR="$1"
OUTPUT_DIR="$2"

[[ -d "$INPUT_DIR" ]] || { echo "Error: input directory not found: $INPUT_DIR" >&2; exit 1; }
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
REPORT="$OUTPUT_DIR/research_report.md"

# Collect the .pdb files to process (top level of INPUT_DIR only).
PDB_FILES=()
for f in "$INPUT_DIR"/*.pdb; do
  [[ -e "$f" ]] && PDB_FILES+=("$f")
done
[[ ${#PDB_FILES[@]} -gt 0 ]] || { echo "Error: no .pdb files found in $INPUT_DIR" >&2; exit 1; }

# The four combinations to run, as "ENGINE REFINEMENT" pairs.
COMBOS=(
  "tinker default"
  "tinker per_layer"
  "gromacs default"
  "gromacs per_layer"
)

# Append a line to the shared report (no timestamps, unlike log()).
report() {
  printf '%s\n' "$*" >> "$REPORT"
}

# ----------------------------------------------------------------------------
# Emit one combo's Markdown block into the report.
#   $1 engine  $2 run_type  $3 work dir  $4 analysis (RUN_DIR) dir
# ----------------------------------------------------------------------------
emit_report_area() {
  local engine="$1"
  local run_type="$2"
  local work="$3"
  local run_dir="$4"
  local log_file="$work/application.LOG"

  report ""
  report "### $engine $run_type"
  report ""

  # --- Total runtime -------------------------------------------------------
  local runtime="N/A"
  if [[ -f "$log_file" ]]; then
    local rt_line
    rt_line="$(grep 'Total runtime was' "$log_file" | tail -1 || true)"
    if [[ -n "$rt_line" ]]; then
      runtime="$(printf '%s' "$rt_line" | grep -oE '[0-9]+ seconds' | tail -1)"
      [[ -n "$runtime" ]] || runtime="N/A"
    fi
  fi
  report "- **Total runtime:** $runtime"

  # --- C-alpha RMSD avg / max ---------------------------------------------
  local rmst="$run_dir/O_system_rmst.txt"
  if [[ -f "$rmst" ]]; then
    local stats
    stats="$(awk 'NR>1 {s+=$2; c++; if(c==1 || $2>m) m=$2} END{if(c) printf "%.3f %.3f", s/c, m}' "$rmst")"
    if [[ -n "$stats" ]]; then
      report "- **C-alpha RMSD:** avg = ${stats%% *}, max = ${stats##* }"
    else
      report "- **C-alpha RMSD:** _no data rows in O_system_rmst.txt_"
    fi
  else
    report "- **C-alpha RMSD:** _O_system_rmst.txt not found_"
  fi

  # --- SR values from the SR-bearing O_system_mt*.lst files ---------------
  # O_system_mtIDc.lst carries no SR value, so it is intentionally excluded.
  # O_system_mtMER.lst is the most important metric and is bold-highlighted.
  local SR_FILES=(
    O_system_mtMER.lst
    O_system_mtIDa.lst
    O_system_mtIDe.lst
    O_system_mtPOS.lst
  )
  report ""
  local found_lst=false
  local name f sr label
  for name in "${SR_FILES[@]}"; do
    f="$run_dir/$name"
    [[ -e "$f" ]] || continue
    if [[ "$found_lst" == false ]]; then
      report "| .lst file | final SR |"
      report "| --- | --- |"
      found_lst=true
    fi
    sr="$(awk '/^SR/{print $2}' "$f" | tail -1)"
    [[ -n "$sr" ]] || sr="N/A"
    if [[ "$name" == "O_system_mtMER.lst" ]]; then
      report "| **$name** | **$sr** |"
    else
      report "| $name | $sr |"
    fi
  done
  if [[ "$found_lst" == false ]]; then
    report "_no SR-bearing O_system_mt*.lst files found_"
  fi
}

# ----------------------------------------------------------------------------
# Report header
# ----------------------------------------------------------------------------
{
  printf '# WDROP research report\n\n'
  printf '_Generated: %(%F %T)T_\n\n' -1
  printf 'Input directory: %s\n\n' "$INPUT_DIR"
  printf 'Output directory: %s\n\n' "$OUTPUT_DIR"
  printf 'PDB files: %s\n' "$(for f in "${PDB_FILES[@]}"; do basename "$f"; done | paste -sd' ' -)"
} > "$REPORT"

log "Research run started. Input dir: $INPUT_DIR"
log "Output dir: $OUTPUT_DIR"
log "Report: $REPORT"
log "PDB files to process: ${#PDB_FILES[@]}"

# ----------------------------------------------------------------------------
# Per-PDB loop
# ----------------------------------------------------------------------------
for SRC_PDB in "${PDB_FILES[@]}"; do
  PDBID="$(basename "${SRC_PDB%.pdb}")"
  PDB_ROOT="$OUTPUT_DIR/$PDBID"
  mkdir -p "$PDB_ROOT"

  log "=== $PDBID ==="
  report ""
  report "## $PDBID"

  # Keep a pristine master copy that is NEVER fed to a run.
  ORIGINAL="$PDB_ROOT/$PDBID.original.pdb"
  cp "$SRC_PDB" "$ORIGINAL"

  for COMBO in "${COMBOS[@]}"; do
    ENGINE="${COMBO%% *}"
    REFINEMENT="${COMBO##* }"
    WORK="$PDB_ROOT/${ENGINE}_${REFINEMENT}"
    mkdir -p "$WORK"

    # Fresh copy of the pristine master for this combo (wdrop.sh mutates its
    # input in place, so nothing can be shared between combos).
    cp "$ORIGINAL" "$WORK/$PDBID.pdb"
    cp "$ORIGINAL" "$WORK/${PDBID}_r.pdb"

    # Analysis dir: default 5 layers -> <PDB>_l5_{REFINEMENT}/ ; one layer is
    # deposited per iteration, so the final iteration lands in the numbered subdir 5.
    RUN_DIR="$WORK/${PDBID}_l5_${REFINEMENT}/5"

    # TINKER_GPU is inherited from the caller's environment (set it manually
    # before running this script to enable Tinker9 GPU mode).
    log "Running $PDBID $ENGINE $REFINEMENT ..."
    if ( cd "$WORK" && "$SCRIPT_DIR/wdrop.sh" \
           "$PDBID.pdb" "$ENGINE" "${PDBID}_r.pdb" --refinement "$REFINEMENT" ); then
      log "DONE: $PDBID $ENGINE $REFINEMENT"
    else
      log "WARNING: $PDBID $ENGINE $REFINEMENT failed (rc=$?); recording partial results."
    fi

    emit_report_area "$ENGINE" "$REFINEMENT" "$WORK" "$RUN_DIR"
  done
done

log "Report written to $REPORT"
