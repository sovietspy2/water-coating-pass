#!/bin/bash
# ============================================================================
# research.sh — batch benchmark + CSV history generator for the WDROP pipeline.
#
# Processes every .pdb file found in an input directory (the files are provided
# and reviewed manually — this script does NOT download anything). For each PDB
# it runs every combo listed in COMBOS below. Each combo is a full wdrop.sh
# invocation: the engine (tinker|gromacs) followed by any wdrop.sh flags, e.g.
#   tinker --layers 5 --intermediate-md-ns 0.1 --final-md-ns 0.5
# Edit COMBOS to define the runs you want.
# Each combo runs in its OWN working directory (the _l{L}_{refinement} output tags
# collide between combos, and each combo needs its own application.LOG for runtime
# extraction). To run Tinker in GPU mode, export TINKER_GPU=1 before invoking
# this script.
#
# Results are appended as rows to a single committed, append-only CSV history
# (testing/research_history.csv by default; override with RESEARCH_CSV). The file
# is NEVER truncated — it is the long-term performance record, meant to be bulk-
# loaded into a DB. One row is written per PDB x combo, with columns:
#   id, pdb_name, engine, layers, refinement, intermediate_md_ns, final_md_ns,
#   total_runtime_seconds, rmsd_avg, rmsd_max, sr_mer, sr_ida, sr_ide, sr_pos,
#   run_timestamp, commit_hash, cpu, gpu
# Metrics come from: total runtime (last "Total runtime was N seconds" line of
# application.LOG), avg/max C-alpha RMSD from O_system_rmst.txt, and the final SR
# value from each O_system_mt{MER,IDa,IDe,POS}.lst. Missing metrics are left empty.
# cpu/gpu come from the CPU/GPU env vars (UNKNOWN when unset).
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
  <OUTPUT_DIR>   directory where per-combo working files are written

Results are appended to \$RESEARCH_CSV (default: testing/research_history.csv),
which is committed and only ever appended to — never replaced or deleted.
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

# Append-only CSV history, committed under testing/. Overridable via RESEARCH_CSV.
# Resolved to an absolute path so a per-combo `cd` into the working dir can't move it.
RESEARCH_CSV="${RESEARCH_CSV:-$SCRIPT_DIR/../testing/research_history.csv}"
mkdir -p "$(dirname "$RESEARCH_CSV")"
RESEARCH_CSV="$(cd "$(dirname "$RESEARCH_CSV")" && pwd)/$(basename "$RESEARCH_CSV")"

CSV_HEADER="id,pdb_name,engine,layers,refinement,intermediate_md_ns,final_md_ns,total_runtime_seconds,rmsd_avg,rmsd_max,sr_mer,sr_ida,sr_ide,sr_pos,run_timestamp,commit_hash,cpu,gpu"

# Collect the .pdb files to process (top level of INPUT_DIR only).
PDB_FILES=()
for f in "$INPUT_DIR"/*.pdb; do
  [[ -e "$f" ]] && PDB_FILES+=("$f")
done
[[ ${#PDB_FILES[@]} -gt 0 ]] || { echo "Error: no .pdb files found in $INPUT_DIR" >&2; exit 1; }

# The combinations to run. Each entry is a full wdrop.sh command line: the engine
# (tinker|gromacs) followed by any wdrop.sh flags (--layers, --intermediate-md-ns,
# --final-md-ns). Flags left out fall back to wdrop.sh defaults (--layers 5,
# --intermediate-md-ns 0 = no intermediate MD, --final-md-ns 0.1). Intermediate MD
# runs iff --intermediate-md-ns > 0. The input/reference PDBs are added
# automatically per combo; do NOT list them here.
COMBOS=(
  "tinker --layers 5 --intermediate-md-ns 0.01 --final-md-ns 0.1"
  "tinker --layers 5 --intermediate-md-ns 0.1 --final-md-ns 0.5"
  "tinker --layers 5 --intermediate-md-ns 0 --final-md-ns 0.1"
  "gromacs --layers 5 --intermediate-md-ns 0.01 --final-md-ns 0.1"
  "gromacs --layers 5 --intermediate-md-ns 0.1 --final-md-ns 0.5"
  "gromacs --layers 5 --intermediate-md-ns 0 --final-md-ns 0.1"
)

# Escape one CSV field (RFC 4180): wrap in double quotes and double any embedded
# quotes, but ONLY when the value contains a comma, quote, or newline. Numbers and
# enum-like tokens pass through untouched; free-text fields (cpu/gpu) are protected.
csv_field() {
  local v="$1"
  case "$v" in
    *[,\"$'\n']*) printf '"%s"' "${v//\"/\"\"}" ;;
    *)            printf '%s' "$v" ;;
  esac
}

# Write the header row once, only when the CSV is missing or empty. Never truncates.
ensure_csv_header() {
  if [[ ! -s "$RESEARCH_CSV" ]]; then
    printf '%s\n' "$CSV_HEADER" >> "$RESEARCH_CSV"
  fi
}

# Next auto-increment id: max existing numeric id in column 1 + 1 (1 when the file
# holds only the header / no data rows). Continues across separate runs.
next_csv_id() {
  local maxid
  maxid="$(awk -F, 'NR>1 && $1 ~ /^[0-9]+$/ {if ($1 > m) m = $1} END {print m + 0}' "$RESEARCH_CSV")"
  printf '%s' "$(( maxid + 1 ))"
}

# Slugify a combo string into a filesystem-safe, unique working-dir name.
# Collapses runs of non-alphanumerics to single '_' and trims leading/trailing '_'
# (e.g. "tinker --layers 5 --intermediate-md-ns 0.1" -> "tinker_layers_5_intermediate_md_ns_0_1").
combo_slug() {
  printf '%s' "$1" | tr -cs '[:alnum:]' '_' | sed 's/^_*//; s/_*$//'
}

# Extract the effective value of a wdrop.sh flag from a combo's arg list, honoring
# both "--flag value" and "--flag=value" forms. Falls back to $2 (the wdrop.sh
# default) when the flag is absent.
#   $1 flag name (e.g. --layers)   $2 default   $3.. the combo args
combo_flag() {
  local flag="$1" def="$2"; shift 2
  local prev=""
  local a
  for a in "$@"; do
    case "$a" in
      "$flag="*) printf '%s' "${a#*=}"; return ;;
    esac
    [[ "$prev" == "$flag" ]] && { printf '%s' "$a"; return; }
    prev="$a"
  done
  printf '%s' "$def"
}

# Final SR value from an O_system_mt*.lst file (last line starting with "SR",
# column 2). Empty string when the file is missing or carries no SR line.
sr_value() {
  local f="$1"
  [[ -f "$f" ]] || { printf ''; return; }
  awk '/^SR/{print $2}' "$f" | tail -1
}

# ----------------------------------------------------------------------------
# Append one combo's result as a CSV row. Uses the run-level constants
# COMMIT_HASH, CPU_INFO, GPU_INFO and the CSV_NEXT_ID counter (incremented here).
#   $1 pdb_name  $2 engine  $3 layers  $4 refinement
#   $5 intermediate_md_ns  $6 final_md_ns  $7 work dir  $8 analysis (RUN_DIR) dir
# ----------------------------------------------------------------------------
append_csv_row() {
  local pdb_name="$1" engine="$2" layers="$3" refinement="$4"
  local intermediate_md_ns="$5" final_md_ns="$6" work="$7" run_dir="$8"
  local log_file="$work/application.LOG"

  # --- Total runtime (seconds) --------------------------------------------
  local runtime=""
  if [[ -f "$log_file" ]]; then
    local rt_line
    rt_line="$(grep 'Total runtime was' "$log_file" | tail -1 || true)"
    if [[ -n "$rt_line" ]]; then
      runtime="$(printf '%s' "$rt_line" | grep -oE '[0-9]+ seconds' | tail -1 | grep -oE '[0-9]+' || true)"
    fi
  fi

  # --- C-alpha RMSD avg / max (empty when no data) ------------------------
  local rmsd_avg="" rmsd_max=""
  local rmst="$run_dir/O_system_rmst.txt"
  if [[ -f "$rmst" ]]; then
    local stats
    stats="$(awk 'NR>1 {s+=$2; c++; if(c==1 || $2>m) m=$2} END{if(c) printf "%.3f %.3f", s/c, m}' "$rmst")"
    if [[ -n "$stats" ]]; then
      rmsd_avg="${stats%% *}"
      rmsd_max="${stats##* }"
    fi
  fi

  # --- Final SR values (O_system_mtIDc.lst carries no SR, so it is excluded) -
  local sr_mer sr_ida sr_ide sr_pos
  sr_mer="$(sr_value "$run_dir/O_system_mtMER.lst")"
  sr_ida="$(sr_value "$run_dir/O_system_mtIDa.lst")"
  sr_ide="$(sr_value "$run_dir/O_system_mtIDe.lst")"
  sr_pos="$(sr_value "$run_dir/O_system_mtPOS.lst")"

  local run_timestamp
  run_timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

  # Build the row through csv_field so free-text fields (cpu/gpu) stay well-formed.
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$(csv_field "$CSV_NEXT_ID")" \
    "$(csv_field "$pdb_name")" \
    "$(csv_field "$engine")" \
    "$(csv_field "$layers")" \
    "$(csv_field "$refinement")" \
    "$(csv_field "$intermediate_md_ns")" \
    "$(csv_field "$final_md_ns")" \
    "$(csv_field "$runtime")" \
    "$(csv_field "$rmsd_avg")" \
    "$(csv_field "$rmsd_max")" \
    "$(csv_field "$sr_mer")" \
    "$(csv_field "$sr_ida")" \
    "$(csv_field "$sr_ide")" \
    "$(csv_field "$sr_pos")" \
    "$(csv_field "$run_timestamp")" \
    "$(csv_field "$COMMIT_HASH")" \
    "$(csv_field "$CPU_INFO")" \
    "$(csv_field "$GPU_INFO")" \
    >> "$RESEARCH_CSV"

  (( CSV_NEXT_ID++ ))
}

# ----------------------------------------------------------------------------
# Run-level constants (captured once, shared by every CSV row)
# ----------------------------------------------------------------------------
COMMIT_HASH="$(git -C "$SCRIPT_DIR" rev-parse HEAD 2>/dev/null || echo UNKNOWN)"
CPU_INFO="${CPU:-UNKNOWN}"
GPU_INFO="${GPU:-UNKNOWN}"

ensure_csv_header
CSV_NEXT_ID="$(next_csv_id)"

log "Research run started. Input dir: $INPUT_DIR"
log "Output dir: $OUTPUT_DIR"
log "CSV history: $RESEARCH_CSV (starting at id $CSV_NEXT_ID)"
log "commit=$COMMIT_HASH cpu=$CPU_INFO gpu=$GPU_INFO"
log "PDB files to process: ${#PDB_FILES[@]}"

# ----------------------------------------------------------------------------
# Per-PDB loop
# ----------------------------------------------------------------------------
for SRC_PDB in "${PDB_FILES[@]}"; do
  PDBID="$(basename "${SRC_PDB%.pdb}")"
  PDB_ROOT="$OUTPUT_DIR/$PDBID"
  mkdir -p "$PDB_ROOT"

  log "=== $PDBID ==="

  # Keep a pristine master copy that is NEVER fed to a run.
  ORIGINAL="$PDB_ROOT/$PDBID.original.pdb"
  cp "$SRC_PDB" "$ORIGINAL"

  for COMBO in "${COMBOS[@]}"; do
    # Split the combo into its tokens: engine + wdrop.sh flags.
    read -ra COMBO_TOKENS <<< "$COMBO"
    [[ ${#COMBO_TOKENS[@]} -ge 1 ]] || { log "Skipping empty combo entry"; continue; }
    ENGINE="${COMBO_TOKENS[0]}"
    WDROP_ARGS=("${COMBO_TOKENS[@]:1}") # everything after the engine

    # Derive the effective params (mirroring wdrop.sh defaults) — used both to
    # locate the output dir and to record them in the CSV row.
    LAYERS="$(combo_flag --layers 5 "${WDROP_ARGS[@]}")"
    INTERMEDIATE_MD_NS="$(combo_flag --intermediate-md-ns 0 "${WDROP_ARGS[@]}")"
    FINAL_MD_NS="$(combo_flag --final-md-ns 0.1 "${WDROP_ARGS[@]}")"
    # REFINEMENT is a derived label (mirrors wdrop.sh): per_layer iff intermediate MD
    # runs (--intermediate-md-ns > 0), else default. Drives the output-dir name below.
    if md_enabled "$INTERMEDIATE_MD_NS"; then REFINEMENT="per_layer"; else REFINEMENT="default"; fi

    # Unique working dir per combo (engine+refinement alone can collide when only
    # the MD lengths differ), so slug the whole combo string.
    SLUG="$(combo_slug "$COMBO")"
    WORK="$PDB_ROOT/$SLUG"
    mkdir -p "$WORK"

    # Fresh copy of the pristine master for this combo (wdrop.sh mutates its
    # input in place, so nothing can be shared between combos).
    cp "$ORIGINAL" "$WORK/$PDBID.pdb"
    cp "$ORIGINAL" "$WORK/${PDBID}_r.pdb"

    # Analysis dir: wdrop.sh writes to <PDB>_l{LAYERS}_{REFINEMENT}/. One layer is
    # deposited per iteration, so the final iteration lands in the numbered subdir
    # {LAYERS} — except a single-layer run, which has no numbered subdir.
    if (( LAYERS > 1 )); then
      RUN_DIR="$WORK/${PDBID}_l${LAYERS}_${REFINEMENT}/${LAYERS}"
    else
      RUN_DIR="$WORK/${PDBID}_l${LAYERS}_${REFINEMENT}"
    fi

    # TINKER_GPU is inherited from the caller's environment (set it manually
    # before running this script to enable Tinker9 GPU mode).
    log "Running $PDBID: $COMBO ..."
    if ( cd "$WORK" && "$SCRIPT_DIR/wdrop.sh" \
           "$PDBID.pdb" "$ENGINE" "${PDBID}_r.pdb" "${WDROP_ARGS[@]}" ); then
      log "DONE: $PDBID $COMBO"
    else
      log "WARNING: $PDBID '$COMBO' failed (rc=$?); recording partial results."
    fi

    append_csv_row "$PDBID" "$ENGINE" "$LAYERS" "$REFINEMENT" \
      "$INTERMEDIATE_MD_NS" "$FINAL_MD_NS" "$WORK" "$RUN_DIR"
  done
done

log "CSV history updated: $RESEARCH_CSV"
