#!/usr/bin/env bash
# Shared plumbing for the minimize/MD determinism suites.
#
# These suites drive the ENGINE BINARIES directly (minimize, dynamic, dynamic9,
# gmx mdrun) rather than pipeline/tinker.sh or pipeline/gromacs.sh. The commands,
# key files and mdp files below are lifted verbatim from those pipeline scripts so
# the measurement reflects what the pipeline actually runs, but invoking the binaries
# directly is what lets us hold every input byte-identical across replicates and
# attribute the variance to one step.
#
# Sourced by:
#   testing/minimize-deterministic-test/{GROMACS,TINKER}/run.sh
#   testing/md-deterministic-test/{GROMACS,TINKER}/run.sh

# ============================================================================
# PATHS
# ============================================================================

DETERMINISM_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly DETERMINISM_LIB_DIR
PROJECT_ROOT="$(cd "$DETERMINISM_LIB_DIR/../.." && pwd)"
readonly PROJECT_ROOT
readonly PIPELINE_DIR="$PROJECT_ROOT/pipeline"
readonly DEFAULT_INPUT_PDB="$PROJECT_ROOT/testing/test_pdbs/1PSV_cryst.pdb"
readonly COMPARE_PY="$DETERMINISM_LIB_DIR/compare_structures.py"
readonly PARSE_SCALARS_PY="$DETERMINISM_LIB_DIR/parse_scalars.py"

# ============================================================================
# COMMON OPTIONS (set by parse_common_args)
# ============================================================================

REPLICATES=5
MD_NS="0.02"
ARMS_SELECTOR=""          # pipeline | hardened | all; per-suite default applied by the caller
KEEP_WORK_DIR=false
INPUT_PDB="$DEFAULT_INPUT_PDB"
TARGET_FRAMES=1000        # matches run_mm_step's hardcoded 1000 in pipeline_common.sh
WORK_DIR=""
RESULTS_DIR=""
RUN_FAILED=false

# ============================================================================
# OUTPUT
# ============================================================================

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

common_usage_options() {
  cat <<'EOF'
  -n, --replicates N     number of identical runs to compare (default: 5)
  -i, --input PATH       input PDB (default: testing/test_pdbs/1PSV_cryst.pdb)
      --md-ns NS         MD length in nanoseconds (MD suite only, default: 0.02)
      --arms WHICH       pipeline | hardened | all
      --keep             keep the temporary work directory after the run
  -h, --help             show this help message
EOF
}

# Consumes the options above and leaves anything it does not recognise in
# DETERMINISM_UNPARSED_ARGS for the caller to handle.
DETERMINISM_UNPARSED_ARGS=()
parse_common_args() {
  DETERMINISM_UNPARSED_ARGS=()

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -n|--replicates)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        REPLICATES="$2"; shift 2 ;;
      -i|--input)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        INPUT_PDB="$2"; shift 2 ;;
      --md-ns)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        MD_NS="$2"; shift 2 ;;
      --arms)
        [[ $# -ge 2 ]] || die "$1 requires a value"
        ARMS_SELECTOR="$2"; shift 2 ;;
      --keep)
        KEEP_WORK_DIR=true; shift ;;
      *)
        DETERMINISM_UNPARSED_ARGS+=("$1"); shift ;;
    esac
  done
}

validate_common_args() {
  [[ "$REPLICATES" =~ ^[2-9][0-9]*$ ]] \
    || die "--replicates must be an integer >= 2 (got '$REPLICATES'); comparing needs at least two runs"
  [[ "$MD_NS" =~ ^[0-9]+(\.[0-9]+)?$ ]] \
    || die "--md-ns must be a non-negative number (got '$MD_NS')"

  case "${ARMS_SELECTOR:-all}" in
    pipeline|hardened|all) ;;
    *) die "--arms must be pipeline, hardened or all (got '$ARMS_SELECTOR')" ;;
  esac
}

# Resolve the arm list for a case from ARMS_SELECTOR.
# $1 = per-suite default selector used when the user passed none.
# $2 = "single" when the case has no hardened counterpart (Tinker GPU).
resolve_arms() {
  local DEFAULT_SELECTOR="$1"
  local SINGLE_ARM_ONLY="${2:-}"
  local SELECTOR="${ARMS_SELECTOR:-$DEFAULT_SELECTOR}"

  if [[ "$SINGLE_ARM_ONLY" == "single" ]]; then
    printf 'pipeline\n'
    return 0
  fi

  case "$SELECTOR" in
    pipeline) printf 'pipeline\n' ;;
    hardened) printf 'hardened\n' ;;
    all)      printf 'pipeline\nhardened\n' ;;
  esac
}

# ============================================================================
# WORK DIRECTORY / INPUT INTEGRITY
# ============================================================================

# The input PDB is a committed fixture and must never be touched. Everything runs
# on a copy inside a temp dir; the source hash is checked again on exit.
INPUT_PDB_SHA_BEFORE=""

require_input_pdb() {
  [[ -f "$INPUT_PDB" ]] || die "input PDB not found: $INPUT_PDB"
  [[ -s "$INPUT_PDB" ]] || die "input PDB is empty (0 bytes): $INPUT_PDB
Restore it before running this suite, e.g. from testing/pdb-preprocess-test/."

  INPUT_PDB="$(cd -P "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"
  INPUT_PDB_SHA_BEFORE="$(sha256sum "$INPUT_PDB" | awk '{print $1}')"
}

verify_input_pdb_untouched() {
  local SHA_AFTER
  SHA_AFTER="$(sha256sum "$INPUT_PDB" 2>/dev/null | awk '{print $1}')" || SHA_AFTER="<missing>"

  if [[ "$SHA_AFTER" != "$INPUT_PDB_SHA_BEFORE" ]]; then
    printf 'ERROR: the input PDB was modified by this run: %s\n' "$INPUT_PDB" >&2
    printf '       before=%s after=%s\n' "$INPUT_PDB_SHA_BEFORE" "$SHA_AFTER" >&2
    return 1
  fi
}

setup_work_dir() {
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/wdrop-determinism.XXXXXX")"
  log "Work directory: $WORK_DIR"
}

# Registered with `trap ... EXIT` by every leaf run.sh.
determinism_cleanup() {
  local STATUS=$?

  if [[ -n "${INPUT_PDB_SHA_BEFORE:-}" ]]; then
    verify_input_pdb_untouched || STATUS=1
  fi

  if [[ -n "${WORK_DIR:-}" && -d "$WORK_DIR" ]]; then
    if [[ "$STATUS" -eq 0 && "$KEEP_WORK_DIR" == false && "$RUN_FAILED" == false ]]; then
      rm -rf -- "$WORK_DIR"
    else
      log "Work directory kept at: $WORK_DIR"
    fi
  fi

  return "$STATUS"
}

# ============================================================================
# PYTHON ENVIRONMENT
# ============================================================================

# The pipeline auto-creates .venv in the repo root; reuse it so pdbfixer is available.
activate_venv() {
  local VENV_DIR="$PROJECT_ROOT/.venv"

  [[ -d "$VENV_DIR" ]] || python3 -m venv "$VENV_DIR"
  [[ -f "$VENV_DIR/bin/activate" ]] || die "missing venv activate script: $VENV_DIR/bin/activate"

  # shellcheck disable=SC1091
  source "$VENV_DIR/bin/activate"
}

# ============================================================================
# FIXTURE
# ============================================================================

# Copy the input PDB into the work dir, run the pipeline's preprocessing, then
# deposit exactly one water layer. wdrop is already proven deterministic
# (testing/algo-deterministic-test), so building this once and sharing it across
# replicates is what keeps the measured variance attributable to the engine.
#
# Echoes the path of the prepared PDB.
build_fixture_pdb() {
  local FIXTURE_DIR="$WORK_DIR/fixture"
  local FIXTURE_PDB="$FIXTURE_DIR/fixture.pdb"
  local WAT_PDB="$FIXTURE_DIR/fixture_1WAT.pdb"

  mkdir -p -- "$FIXTURE_DIR"
  cp -- "$INPUT_PDB" "$FIXTURE_PDB"

  log "Fixture step 1/3: pdb-preprocessor.py --target (single model, repair atoms, strip waters)"
  "$PIPELINE_DIR/pdb-preprocessor.py" --target "$FIXTURE_PDB" >"$FIXTURE_DIR/preprocess.log" 2>&1 \
    || { cat "$FIXTURE_DIR/preprocess.log" >&2; die "pdb-preprocessor.py failed"; }

  log "Fixture step 2/3: wdrop --layers 1 (deposit one water layer)"
  ( cd "$FIXTURE_DIR" \
      && wdrop --file fixture.pdb --sigma 1.8 --weed-dist 3.5 --layers 1 ) \
    >"$FIXTURE_DIR/wdrop.log" 2>&1 \
    || { cat "$FIXTURE_DIR/wdrop.log" >&2; die "wdrop failed"; }
  [[ -f "$WAT_PDB" ]] || die "wdrop did not produce the expected output: $WAT_PDB"

  log "Fixture step 3/3: format_pdb.py (canonical columns, OW->O, SOL->WAT)"
  "$PIPELINE_DIR/format_pdb.py" "$WAT_PDB" >"$FIXTURE_DIR/format.log" 2>&1 \
    || { cat "$FIXTURE_DIR/format.log" >&2; die "format_pdb.py failed"; }

  log "Fixture ready: $WAT_PDB ($(grep -c '^ATOM\|^HETATM' "$WAT_PDB") atoms)"
  printf '%s\n' "$WAT_PDB"
}

# ============================================================================
# UNIT CONVERSION (lifted from pipeline/tinker.sh:7-34 and pipeline/gromacs.sh:4-28)
# ============================================================================

# The locals here are underscore-prefixed on purpose: callers legitimately hold
# readonly DT_FS / DT_PS constants of their own, and a plain `local DT_FS` would
# abort with "readonly variable".
ns_to_steps() {
  local _NS="${1:?usage: ns_to_steps <nanoseconds> <dt_fs>}"
  local _DT_FS="${2:-1.0}"
  awk -v NS="$_NS" -v DT="$_DT_FS" 'BEGIN { printf "%.0f\n", (NS * 1000000) / DT }'
}

ns_to_ps() {
  local _NS="${1:?usage: ns_to_ps <nanoseconds>}"
  awk -v NS="$_NS" 'BEGIN { printf "%.6f\n", NS * 1000.0 }'
}

steps_to_ps() {
  local _STEPS="${1:?usage: steps_to_ps <steps> <dt_ps>}"
  local _DT_PS="${2:-0.001}"
  awk -v STEPS="$_STEPS" -v DT="$_DT_PS" 'BEGIN { printf "%.6f\n", STEPS * DT }'
}

# Frames are saved every N_STEPS/TARGET_FRAMES steps, floored at 1.
save_every_steps() {
  local _N_STEPS="$1"
  local _FRAMES="$2"
  local _EVERY=$(( _N_STEPS / _FRAMES ))
  (( _EVERY < 1 )) && _EVERY=1
  printf '%d\n' "$_EVERY"
}

# ============================================================================
# TINKER INPUTS (lifted from pipeline/tinker.sh)
# ============================================================================

# tinker.sh:113 — protein atoms are every XYZ atom before the first water (OW/HW).
tinker_protein_atom_count() {
  local XYZ_FILE="$1"
  awk 'NR>1 && ($2=="OW" || $2=="HW") {print NR-2; exit}' "$XYZ_FILE"
}

# tinker.sh:103-108 — PDB -> XYZ. Runs in the given directory; echoes the .xyz path.
tinker_build_xyz() {
  local DIR="$1"
  local PDB_PATH="$2"
  local BASE
  BASE="$(basename "${PDB_PATH%.pdb}")"

  cp -f -- "$PIPELINE_DIR/amber99.prm" "$DIR/amber99.prm"
  cp -f -- "$PDB_PATH" "$DIR/$BASE.pdb"

  ( cd "$DIR" && pdbxyz "$BASE.pdb" >pdbxyz.log 2>&1 <<EOF
ALL
amber99.prm
EOF
  ) || { cat "$DIR/pdbxyz.log" >&2; die "pdbxyz failed"; }

  [[ -f "$DIR/$BASE.xyz" ]] || die "pdbxyz did not produce $DIR/$BASE.xyz"
  printf '%s\n' "$DIR/$BASE.xyz"
}

# tinker.sh:129-133. $2 = protein atom count, $3 = arm.
# The pipeline ships OPENMP-THREADS 1 here; the `hardened` arm drops it to expose the
# thread-scheduling non-determinism that line was added to suppress.
tinker_write_minimize_key() {
  local KEY_PATH="$1"
  local PROTEIN_ATOMS="$2"
  local ARM="$3"

  cat > "$KEY_PATH" <<EOF
RESTRAIN-POSITION -1 ${PROTEIN_ATOMS} 2.0
PARAMETERS amber99.prm
EOF

  if [[ "$ARM" == "pipeline" ]]; then
    printf 'OPENMP-THREADS 1\n' >> "$KEY_PATH"
  fi
}

# tinker.sh:152-171. $3 = arm, $4 = "gpu" to add the Tinker9 lines.
# The pipeline ships md.key WITHOUT a thread pin, so here `hardened` ADDS
# OPENMP-THREADS 1 — the one-line change TINKER-DETERMINISM.md section 4 leaves open.
tinker_write_md_key() {
  local KEY_PATH="$1"
  local PROTEIN_ATOMS="$2"
  local ARM="$3"
  local GPU_MODE="${4:-}"

  cat > "$KEY_PATH" <<EOF
PARAMETERS amber99.prm
RANDOMSEED 28480426
RESTRAIN-POSITION -1 ${PROTEIN_ATOMS} 300.0

RATTLE
vdw-cutoff 9.0
chg-cutoff 9.0
ARCHIVE
EOF

  if [[ "$GPU_MODE" == "gpu" ]]; then
    printf 'INTEGRATOR VERLET\nREMOVE-INERTIA 0\n' >> "$KEY_PATH"
  fi

  if [[ "$ARM" == "hardened" ]]; then
    printf 'OPENMP-THREADS 1\n' >> "$KEY_PATH"
  fi
}

# ============================================================================
# GROMACS INPUTS (mdp heredocs copied verbatim from pipeline/gromacs.sh:75-171)
# ============================================================================

gmx_write_st_mdp() {
  local OUTFILE="$1"
  local NSTEPS="${2:-50000}"

  cat > "$OUTFILE" <<EOF
; GROMACS Minimization in Vacuum (Pseudo-PBC)
;
define              = -DPOSRES -DPOSRES_WATER ; Use if you have position restraints
integrator          = steep                    ; Steepest descent minimization
nsteps              = ${NSTEPS}                    ; Max steps (use more than 50k for safety)

; Neighbor searching parameters for pseudo-PBC
nstlist             = 10
cutoff-scheme       = Verlet
ns_type             = simple                   ; Simpler neighbor search for EM
rlist               = 333.3                    ; Neighbor list cutoff (must be large)

; Electrostatics and VdW for vacuum
coulombtype         = cutoff                   ; Simple cutoff
rcoulomb            = 333.3                    ; Large cutoff simulates vacuum (must be < 0.5 * Box Vector)
vdwtype             = cutoff                   ; Simple cutoff
rvdw                = 333.3                    ; Large cutoff simulates vacuum

; Energy minimizing settings
emtol               = 200.0                    ; Stop when the maximum force is below 200 kJ/mol/nm
emstep              = 0.01                     ; Initial step-size
EOF
}

gmx_write_cg_mdp() {
  local OUTFILE="$1"
  local NSTEPS="${2:-250000}"

  cat > "$OUTFILE" <<EOF
; Modified gromacs-cg.mdp for Pseudo-PBC
define              = -DPOSRES -DPOSRES_WATER
constraints         = none
integrator          = cg          ; Conjugate Gradient minimization
nsteps              = ${NSTEPS}      ; Max steps
nstlist             = 10
cutoff-scheme       = Verlet
ns_type             = simple      ; Use simple for minimization
rlist               = 333.3       ; Large cutoff for vacuum
coulombtype         = cutoff      ; Use cutoff with large radius for vacuum
rcoulomb            = 333.3       ; Large cutoff for vacuum
vdwtype             = cutoff      ; Use cutoff with large radius for vacuum
rvdw                = 333.3       ; Large cutoff for vacuum
emtol               = 750         ; Max force target
emstep              = 0.01        ;
EOF
}

gmx_write_md_mdp() {
  local OUTFILE="$1"
  local NSTEPS="${2:-500000}"
  local SAVE_EVERY_STEPS="${3:-500}"
  local DT="${4:-0.002}"

  cat > "$OUTFILE" <<EOF
; Modified gromacs-md.mdp for Pseudo-PBC
cpp                 = /usr/bin/cpp
define              = -DPOSRES
constraints         = all-bonds
integrator          = md
dt                  = ${DT}
nsteps              = ${NSTEPS}
nstcomm             = 500
nstxout             = ${SAVE_EVERY_STEPS}
nstvout             = ${SAVE_EVERY_STEPS}
nstfout             = 0
nstlog              = 500
nstenergy           = 500
nstlist             = 10
ns_type             = grid
cutoff-scheme       = Verlet

; Pseudo-PBC / Vacuum Settings
coulombtype         = cutoff        ; Use cutoff for vacuum
rlist               = 333.3         ; Large cutoff for vacuum
rcoulomb            = 333.3         ; Large cutoff for vacuum
rvdw                = 333.3         ; Large cutoff for vacuum

; Temperature Coupling
Tcoupl              = v-rescale
tc-grps             = Protein non-Protein
tau_t               = 0.1 0.1
ref_t               = 300 300
energygrps          = Protein non-Protein

; NO Pressure Coupling for vacuum simulation

gen_vel             = yes           ; Generate initial velocities
gen_temp            = 300.0
gen_seed            = 28480426       ; Fixed seed
ld_seed             = 28480426
EOF
}

# gromacs.sh:191-195 — topology + pseudo-PBC box. Runs in DIR; leaves conf_largebox.gro,
# topol.top and posre.itp behind for grompp.
gmx_build_topology() {
  local DIR="$1"
  local PDB_PATH="$2"

  cp -f -- "$PDB_PATH" "$DIR/input.pdb"

  ( cd "$DIR" \
      && gmx pdb2gmx -water tip3p -ff amber99sb-ildn -ignh -f input.pdb -o conf.gro \
      && gmx editconf -f conf.gro -o conf_largebox.gro -c -d 0 -bt cubic -box 1000 \
  ) >"$DIR/topology.log" 2>&1 \
    || { tail -40 "$DIR/topology.log" >&2; die "gmx pdb2gmx/editconf failed (see $DIR/topology.log)"; }

  [[ -f "$DIR/conf_largebox.gro" ]] || die "gmx editconf did not produce $DIR/conf_largebox.gro"
}

# Write a per-atom "<resname> <atomname>" file from a .gro.
#
# `gmx trjconv -o x.g96` emits POSITIONRED blocks -- bare coordinates with no atom or
# residue names -- so the protein/water split cannot come from the compared file and
# has to be taken from the topology instead. pdb2gmx's conf.gro has the same atom
# count and order as the .tpr and every trajectory derived from it.
gmx_write_atom_labels() {
  local GRO="$1"
  local OUT="$2"

  awk '
    NR == 2 { N = $1 + 0; next }
    NR > 2 && NR <= N + 2 {
      RES = substr($0, 6, 5); ATOM = substr($0, 11, 5)
      gsub(/^ +| +$/, "", RES); gsub(/^ +| +$/, "", ATOM)
      print RES, ATOM
    }
  ' "$GRO" > "$OUT"

  [[ -s "$OUT" ]] || die "could not derive atom labels from $GRO"
}

# The `hardened` arm for GROMACS. -reprod disables all run-to-run-varying performance
# optimisations; pinning to a single rank/thread removes the reduction-order variance
# that is the Tinker root cause, so the two engines are probed the same way.
gmx_mdrun_arm_flags() {
  local ARM="$1"
  if [[ "$ARM" == "hardened" ]]; then
    printf -- '-reprod\n-ntmpi\n1\n-ntomp\n1\n'
  fi
}

# ============================================================================
# RESULTS
# ============================================================================

# $1 = results root (suite dir), $2 = case slug. Echoes the created directory.
make_results_dir() {
  local RESULTS_ROOT="$1"
  local CASE_SLUG="$2"
  local STAMP
  STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

  local DIR="$RESULTS_ROOT/${STAMP}_${CASE_SLUG}"
  mkdir -p -- "$DIR"
  printf '%s\n' "$DIR"
}

# Command prefix that records CPU utilisation for the command under test.
#
# This exists because a "BIT-IDENTICAL" verdict has two very different causes, and the
# report must let you tell them apart: the engine may be genuinely reproducible, or it
# may simply never have left a single thread — in which case there was no reduction
# order to vary and the verdict says nothing about larger systems. Measured here:
# Tinker `dynamic` stays near 1.0 cores even at 1644 atoms, while `minimize` does not.
#
# Expands to nothing when /usr/bin/time is unavailable; callers must tolerate that.
timing_prefix() {
  local OUT="$1"
  if [[ -x /usr/bin/time ]]; then
    printf '/usr/bin/time\n-f\n%%P %%e\n-o\n%s\n' "$OUT"
  fi
}

# Mean CPU cores used across a set of `timing_prefix` output files ("NNN% SS.SS").
# Echoes nothing when there is no usable timing data.
mean_cpu_cores() {
  [[ $# -gt 0 ]] || return 0
  cat -- "$@" 2>/dev/null \
    | awk '/%/ { gsub(/%/, "", $1); TOTAL += $1; N++ } END { if (N) printf "%.2f\n", TOTAL / (100 * N) }'
}

# Append the measured CPU utilisation for one arm to a report.
# $1 = report path, $2 = label, $3 = command name, $4.. = files from `timing_prefix`.
#
# No threshold is applied on purpose. Tinker `minimize` was measured diverging at only
# 1.09 average cores, so "barely above one core" does NOT imply determinism -- brief
# parallel regions are enough to reorder a reduction. The number is reported as
# context for the verdict, not as a predictor of it.
report_cpu_utilisation() {
  local REPORT_PATH="$1"
  local LABEL="$2"
  local COMMAND_NAME="$3"
  shift 3

  local CORES
  CORES="$(mean_cpu_cores "$@")"
  [[ -n "$CORES" ]] || return 0

  log "CPU utilisation: $COMMAND_NAME averaged $CORES cores (of $(nproc) available)"
  {
    printf '**CPU utilisation — %s: `%s` averaged %s cores of %s available.**\n\n' \
      "$LABEL" "$COMMAND_NAME" "$CORES" "$(nproc)"
    printf '> Thread count is the mechanism behind Tinker non-determinism: OpenMP reduction\n'
    printf '> order can only vary when more than one thread participates. A near-1.0 average\n'
    printf '> does not by itself prove determinism, but a BIT-IDENTICAL verdict at that level\n'
    printf '> says little about a larger system that engages more threads.\n\n'
  } >> "$REPORT_PATH"
}

# Start a report. Extra "key: value" lines are appended to the configuration block.
write_report_header() {
  local REPORT="$1"
  local TITLE="$2"
  shift 2

  {
    printf '# %s\n\n' "$TITLE"
    printf -- '- Generated: %s\n' "$(date -u +'%Y-%m-%d %H:%M:%S UTC')"
    printf -- '- Machine: %s\n' "$(machine_description)"
    printf -- '- Commit: %s\n' "$(git -C "$PROJECT_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    printf -- '- Input PDB: `%s` (sha256 %s)\n' "$INPUT_PDB" "${INPUT_PDB_SHA_BEFORE:0:16}"
    printf -- '- Replicates per arm: %s\n' "$REPLICATES"
    local LINE
    for LINE in "$@"; do
      printf -- '- %s\n' "$LINE"
    done
    printf '\n'
  } > "$REPORT"
}

# Record a one-line verdict. When the suite driver exports DETERMINISM_SUMMARY_FILE
# the line is also collected there for the cross-case summary.
emit_headline() {
  local TEXT="$1"
  log "RESULT: $TEXT"
  if [[ -n "${DETERMINISM_SUMMARY_FILE:-}" ]]; then
    printf '%s\n' "$TEXT" >> "$DETERMINISM_SUMMARY_FILE"
  fi
}

# One line describing the machine, so reports from different boxes are never confused.
# TINKER-DETERMINISM.md section 1: Tinker reproducibility "will not hold across
# different machine types", so results are only comparable within one host.
machine_description() {
  local CPU GPU
  CPU="$(lscpu 2>/dev/null | awk -F: '/^Model name/{gsub(/^ +/,"",$2);print $2;exit}')"
  CPU="${CPU:-unknown}"
  if command -v nvidia-smi >/dev/null 2>&1; then
    GPU="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1)"
  fi
  GPU="${GPU:-none}"

  printf 'host=%s cpu=%s threads=%s gpu=%s\n' "$(hostname)" "$CPU" "$(nproc)" "$GPU"
}
