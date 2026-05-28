#!/bin/bash
set -euo pipefail

ns_to_nsteps() {
  local ns="${1:?usage: ns_to_nsteps <ns>}"
  local dt_ps="${2:-0.002}"

  awk -v ns="$ns" -v dt="$dt_ps" 'BEGIN {
    printf "%.0f", (ns * 1000) / dt
  }'
}

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

command -v gmx >/dev/null || { echo "Error: gmx not found" >&2; exit 1; }
command -v python3 >/dev/null || { echo "Error: python not found" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

INPUT_PDB="$1"
INPUT_ABS="$(cd -P "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"
INPUT_DIR="$(dirname "$INPUT_ABS")"
MD_DURATION="$2"
MOBYWAT_OUTPUT_ENABLED="$3"

# Setup logging in the INPUT_DIR (output directory)
LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging $LOGFILE

log "Starting gromacs.sh (GROMACS mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_ABS=$INPUT_ABS"
log "INPUT_DIR=$INPUT_DIR"
log "MD_DURATION=$MD_DURATION ns"
log "MOBYWAT_OUTPUT_ENABLED=$MOBYWAT_OUTPUT_ENABLED"

touch "${INPUT_DIR}/GROMACS.protocol"
log "Created GROMACS.protocol file"

# --- TIMER SETUP ---
SECONDS=0
TIMEFILE="${INPUT_ABS}-mm-process-time.txt"

## GROMACS parameter helpers

write_cg_mdp() {
  local outfile="$INPUT_DIR/gromacs-cg.mdp"
  local nsteps="${1:-250000}"

  cat > "$outfile" <<EOF
; Modified gromacs-cg.mdp for Pseudo-PBC
define              = -DPOSRES -DPOSRES_WATER
constraints         = none
integrator          = cg          ; Conjugate Gradient minimization
nsteps              = ${nsteps}      ; Max steps
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
  local outfile="$INPUT_DIR/gromacs-md.mdp"
  local nsteps="${1:-500000}" # default is total 1 ns.
  local dt="${2:-0.002}"

  cat > "$outfile" <<EOF
; Modified gromacs-md.mdp for Pseudo-PBC
cpp                 = /usr/bin/cpp
define              = -DPOSRES
constraints         = all-bonds
integrator          = md
dt                  = ${dt}
nsteps              = ${nsteps}
nstcomm             = 500
nstxout             = 500
nstvout             = 500
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
  local outfile="$INPUT_DIR/gromacs-st.mdp"
  local nsteps="${1:-50000}"

  cat > "$outfile" <<EOF
; GROMACS Minimization in Vacuum (Pseudo-PBC)
;
define              = -DPOSRES -DPOSRES_WATER ; Use if you have position restraints
integrator          = steep                    ; Steepest descent minimization
nsteps              = ${nsteps}                    ; Max steps (use more than 50k for safety)

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
  local elapsed="$SECONDS"
  printf 'runtime for mm is %s seconds\n' "$elapsed" > "$TIMEFILE"
  log "Runtime written to $TIMEFILE: $elapsed seconds"
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
log "Step 4: Running grompp for conjugate gradient (using gromacs-cg.mdp)"
run_step gmx grompp -v -f "$INPUT_DIR/gromacs-cg.mdp" -c after_em.gro -r after_em.gro -o cg -p topol.top

log "Step 4: Running mdrun for conjugate gradient"
run_step gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log

# --- 5. MOLECULAR DYNAMICS (MD) ---
STEPS="$(ns_to_nsteps $MD_DURATION)"
log "Creating parameter file: gromacs-md.mdp"
write_md_mdp $STEPS
log "Step 5: Running grompp for molecular dynamics (using gromacs-md.mdp)"
run_step gmx grompp -f "$INPUT_DIR/gromacs-md.mdp" -o md -c after_cg.gro -r after_cg.gro -p topol.top -maxwarn 1

log "Step 5: Running mdrun for molecular dynamics"
run_step gmx mdrun -v -s md -o md.trr -c after_md.gro -g md.log

# --- 6. POST-PROCESSING (Cleaning up the trajectory) ---
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
END_PS="$(awk -v ns="$MD_DURATION" 'BEGIN { printf "%.3f", ns * 1000 }')" # We need this because we introduced time params
run_step gmx trjconv -f system_compact.xtc -s md.tpr -o lastframe_drop.pdb -b "$END_PS" -e "$END_PS" <<EOF
0
EOF

# input: lastframe_drop.pdb output: filtered_renum.pdb
# Removing waters too far away from protein
log "Step 7: Removing extra waters not connected to protein"
run_step python3 "$SCRIPT_DIR/clean_pdb.py" "$INPUT_DIR/lastframe_drop.pdb"

log "Step 8: Post processing PDB"
run_step "$SCRIPT_DIR"/format-pdb.sh "$INPUT_DIR/filtered_renum.pdb"

log "Step 9: Post processing, creating MobyWat compatible file trajectory file."
cp system_compact.xtc mobywat_input.xtc

log "gromacs.sh completed successfully"
