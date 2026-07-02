#!/bin/bash
set -euo pipefail

ns_to_nsteps() {
  local NS="${1:?usage: ns_to_nsteps <ns>}"
  local DT_PS="${2:-0.002}"

  awk -v ns="$NS" -v dt="$DT_PS" 'BEGIN {
    printf "%.0f", (ns * 1000) / dt
  }'
}

ns_to_ps() {
  local NS="${1:?usage: ns_to_ps <ns>}"

  awk -v ns="$NS" 'BEGIN {
    printf "%.6f", ns * 1000.0
  }'
}

steps_to_ps() {
  local STEPS="${1:?usage: steps_to_ps <steps>}"
  local DT_PS="${2:-0.002}"

  awk -v steps="$STEPS" -v dt="$DT_PS" 'BEGIN {
    printf "%.6f", steps * dt
  }'
}

if [[ $# -lt 3 || $# -gt 5 ]]; then
  echo "Usage: $0 <input.pdb> <md_duration_ns> <mobywat_output_enabled> <target_frames> <reference.pdb>" >&2
  exit 1
fi

command -v gmx >/dev/null || { echo "Error: gmx not found" >&2; exit 1; }
command -v python3 >/dev/null || { echo "Error: python not found" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

INPUT_PDB="$1"
INPUT_ABS="$(cd -P "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"
INPUT_DIR="$(dirname "$INPUT_ABS")"
PDB_NAME="$(basename "${INPUT_ABS%.pdb}")"
MD_DURATION="$2"
MOBYWAT_OUTPUT_ENABLED="$3"
TARGET_FRAMES="${4:-1000}"
REFERENCE_PDB="${5:-}"
DT_PS="0.001" # 1 femtosecond

# Setup logging in the INPUT_DIR (output directory)
LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging "$LOGFILE"

log "Starting gromacs.sh (GROMACS mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_ABS=$INPUT_ABS"
log "INPUT_DIR=$INPUT_DIR"
log "MD_DURATION=$MD_DURATION ns"
log "MOBYWAT_OUTPUT_ENABLED=$MOBYWAT_OUTPUT_ENABLED"
log "TARGET_FRAMES=$TARGET_FRAMES"
log "DT_PS=$DT_PS"

touch "${INPUT_DIR}/GROMACS.protocol"
log "Created GROMACS.protocol file"

# --- TIMER SETUP ---
SECONDS=0
TIMEFILE="${INPUT_DIR}/${PDB_NAME}-mm-process-time.txt"

## GROMACS parameter helpers

write_cg_mdp() {
  local OUTFILE="$INPUT_DIR/gromacs-cg.mdp"
  local NSTEPS="${1:-250000}"

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

write_md_mdp() {
  local OUTFILE="$INPUT_DIR/gromacs-md.mdp"
  local NSTEPS="${1:-500000}" # default is total 1 ns.
  local SAVE_EVERY_STEPS="${2:-500}"
  local DT="${3:-0.002}"

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
gen_seed            = 28480426
EOF
}

write_st_mdp() {
  local OUTFILE="$INPUT_DIR/gromacs-st.mdp"
  local NSTEPS="${1:-50000}"

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


##

write_runtime() {
  local ELAPSED="$SECONDS"
  printf 'runtime for mm is %s seconds\n' "$ELAPSED" > "$TIMEFILE"
  log "Runtime written to $TIMEFILE: $ELAPSED seconds"
}
trap write_runtime EXIT

cd "$INPUT_DIR"
log "Changed working directory to: $INPUT_DIR"

log "Step 0: Reformatting PDB"
run_step "$SCRIPT_DIR"/format-pdb.sh "$INPUT_PDB"

# --- 1. SETUP ---
log "Step 1: Running pdb2gmx (generate topology and initial structure)"
run_step gmx pdb2gmx -water tip3p -ff amber99sb-ildn -ignh -f "$INPUT_ABS" -o conf.gro

# --- 2. PREPARATION FOR PSEUDO-PBC ---
log "Step 2: Running editconf (prepare large cubic box for pseudo-PBC)"
run_step gmx editconf -f conf.gro -o conf_largebox.gro -c -d 0 -bt cubic -box 1000

# --- 3. ENERGY MINIMIZATION (EM) ---
log "Creating parameter file: gromacs-st.mdp"
write_st_mdp
log "Step 3: Running grompp for energy minimization (using gromacs-st.mdp)"
run_step gmx grompp -v -f "$INPUT_DIR/gromacs-st.mdp" -c conf_largebox.gro -r conf_largebox.gro -o em -p topol.top -maxwarn 2

log "Step 3: Running mdrun for energy minimization"
run_step gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log

# --- 4. CONSTRAINED MINIMIZATION / CONJUGATE GRADIENT (CG) ---
log "Creating parameter file: gromacs-cg.mdp"
write_cg_mdp 250000 # number of steps
log "Step 4a: Running grompp for conjugate gradient (using gromacs-cg.mdp)"
run_step gmx grompp -v -f "$INPUT_DIR/gromacs-cg.mdp" -c after_em.gro -r after_em.gro -o cg -p topol.top

log "Step 4b: Running mdrun for conjugate gradient"
run_step gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log

# --- 5. MOLECULAR DYNAMICS ---
N_STEPS="$(ns_to_nsteps "$MD_DURATION" "$DT_PS")"
TOTAL_TIME_PS="$(ns_to_ps "$MD_DURATION")"

if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
    SAVE_EVERY_STEPS=$(( N_STEPS / TARGET_FRAMES ))
    if (( SAVE_EVERY_STEPS < 1 )); then
        SAVE_EVERY_STEPS=1
    fi

    ACTUAL_FRAMES=$(( (N_STEPS + SAVE_EVERY_STEPS - 1) / SAVE_EVERY_STEPS ))
    SAVE_INTERVAL_PS="$(steps_to_ps "$SAVE_EVERY_STEPS" "$DT_PS")"
else
    SAVE_EVERY_STEPS=0
    ACTUAL_FRAMES=0
    SAVE_INTERVAL_PS="0.000000"
fi

log "Step 5: Running molecular dynamics"
log "MD_DURATION=$MD_DURATION ns"
log "TOTAL_TIME_PS=$TOTAL_TIME_PS ps"
log "N_STEPS=$N_STEPS"
log "SAVE_EVERY_STEPS=$SAVE_EVERY_STEPS steps"
log "SAVE_INTERVAL_PS=$SAVE_INTERVAL_PS ps"
log "Expected saved frames ~= $ACTUAL_FRAMES"
log "Creating parameter file: gromacs-md.mdp"
write_md_mdp "$N_STEPS" "$SAVE_EVERY_STEPS" "$DT_PS"
log "Step 5a: Running grompp for molecular dynamics (using gromacs-md.mdp)"
run_step gmx grompp -f "$INPUT_DIR/gromacs-md.mdp" -o md -c after_cg.gro -r after_cg.gro -p topol.top -maxwarn 1

log "Step 5b: Running mdrun for molecular dynamics"
if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
  run_step gmx mdrun -v -s md -o md.trr -c after_md.gro -g md.log
else
  run_step gmx mdrun -v -s md -c after_md.gro -g md.log
fi

# --- 6. POST-PROCESSING (Cleaning up the trajectory) ---
if [[ "$MOBYWAT_OUTPUT_ENABLED" == "true" ]]; then
  log "Step 6a: Running trjconv for PBC correction (pbc whole)"
  run_step gmx trjconv -f md.trr -s md.tpr -o pbc_whole.xtc -pbc whole <<EOF
0
EOF

  log "Step 6b: Running trjconv for compact system output"
  run_step gmx trjconv -f pbc_whole.xtc -s md.tpr -o system_compact.xtc -center -pbc mol -ur compact -boxcenter zero <<EOF
1
0
EOF

  log "Step 6c: Creating final frame PDB file"
  END_PS="$TOTAL_TIME_PS" # We need this because we introduced time params
  run_step gmx trjconv -f system_compact.xtc -s md.tpr -o lastframe_drop.pdb -b "$END_PS" -e "$END_PS" <<EOF
0
EOF
else
  log "Step 6: Creating final frame PDB file"
  run_step gmx editconf -f after_md.gro -o lastframe_drop.pdb
fi

# input: lastframe_drop.pdb output: next_step.pdb
# Removing waters too far away from protein
log "Step 7: Removing extra waters not connected to protein"
run_step python3 "$SCRIPT_DIR/remove_unbound_water.py" "$INPUT_DIR/lastframe_drop.pdb" "next_step.pdb"

log "Step 8: Post processing PDB"
run_step "$SCRIPT_DIR"/format-pdb.sh "$INPUT_DIR/next_step.pdb"

if [[ "$MOBYWAT_OUTPUT_ENABLED" != "true" ]]; then
  log "MobyWat output disabled; skipping trajectory PDB generation"
  log "gromacs.sh completed successfully"
  exit 0
fi

#log "Step 9: Post processing, creating MobyWat compatible file trajectory file."
#cp system_compact.xtc mobywat_input.xtc
#cp after_md.gro mobywat_input.gro


log "Step 9: Running mobywat"
if [[ -n "${REFERENCE_PDB:-}" ]]; then

  SYSTEM_REF_PDB="system_ref.pdb"
  cp ${REFERENCE_PDB} ${SYSTEM_REF_PDB}

  log "Making sure Reference PDB is compatible with mobywat."
  run_step "${SCRIPT_DIR}/pdb-preprocessor.py" --reference ${SYSTEM_REF_PDB}


  run_step gmx trjconv -f md.trr -s md.tpr -o pbc1.xtc -pbc whole << EOF
0
EOF

  run_step gmx trjconv -f pbc1.xtc -s md.tpr -o pbc2.xtc -pbc cluster << EOF
1
0
EOF

  run_step gmx trjconv -f pbc2.xtc -s md.tpr -o pbc3.xtc -center -pbc mol -ur compact << EOF
1
0
EOF

  run_step gmx confrms -one -f1 system_ref.pdb -f2 md.tpr -o fit.pdb << EOF
3
3
EOF

  run_step gmx editconf -label A -f fit.pdb -o fit.pdb

  run_step gmx trjconv -f pbc3.xtc -s fit.pdb -o system.xtc -fit progressive << EOF
3
0
EOF

  run_step gmx trjconv -f pbc3.xtc -s fit.pdb -o system_tpy.pdb -b 0 -e 0 -fit progressive << EOF
3
0
EOF

  log "REFERENCE_PDB is present and non-empty: $REFERENCE_PDB, VALIDATION MODE!"
  run_step "${SCRIPT_DIR}"/add-mobywat-analysis-params.sh ${SYSTEM_REF_PDB}

  log "Remove TER operator ID if present from ${SYSTEM_REF_PDB}"
  run_step "${SCRIPT_DIR}"/remove-ter-id.sh ${SYSTEM_REF_PDB}

  log "Running mobywat validation"

else
  log "REFERENCE_PDB is missing or empty, PREDICTION MODE!"
  log "Running mobywat prediction"
  run_step gmx trjconv -f md.trr -s md.tpr -o pbc1.xtc -pbc whole << EOF
0
EOF

  run_step gmx trjconv -f pbc1.xtc -s md.tpr -o pbc2.xtc -pbc cluster << EOF
1
0
EOF

  run_step gmx trjconv -f pbc2.xtc -s md.tpr -o pbc3.xtc -center -pbc mol -ur compact << EOF
1
0
EOF

  run_step gmx trjconv -f pbc3.xtc -s md.tpr -o system.xtc -fit progressive << EOF
3
0
EOF

  run_step gmx trjconv -f pbc3.xtc -s md.tpr -o system_tpy.pdb -b 0 -e 0 -fit progressive << EOF
3
0
EOF

  run_step gmx editconf -label A -f system_tpy.pdb -o system_tpy.pdb # not sure if this is the right way to do it
fi

run_step mobywat -t [A] -w Auto -n 0-1000 -m Prediction -cls MER -v Diagnostic

log "gromacs.sh completed successfully"