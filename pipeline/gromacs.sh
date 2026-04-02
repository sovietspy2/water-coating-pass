#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

command -v gmx >/dev/null || { echo "Error: gmx not found" >&2; exit 1; }
command -v python >/dev/null || { echo "Error: python not found" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline_common.sh"

INPUT_PDB="$1"
INPUT_ABS="$(cd -P "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"
INPUT_DIR="$(dirname "$INPUT_ABS")"

# Setup logging in the INPUT_DIR (output directory)
LOGFILE="${INPUT_DIR}/application.LOG"
setup_logging $LOGFILE

log "Starting gromacs.sh (GROMACS mode)"
log "INPUT_PDB=$INPUT_PDB"
log "INPUT_ABS=$INPUT_ABS"
log "INPUT_DIR=$INPUT_DIR"

touch "${INPUT_DIR}/GROMACS.protocol"
log "Created GROMACS.protocol file"

# --- TIMER SETUP ---
SECONDS=0
TIMEFILE="${INPUT_ABS}-mm-process-time.txt"

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
log "Step 3: Running grompp for energy minimization (using gromacs-st.mdp)"
run_step gmx grompp -v -f "$SCRIPT_DIR/gromacs-st.mdp" -c conf_largebox.gro -r conf_largebox.gro -o em -p topol.top -maxwarn 2

log "Step 3: Running mdrun for energy minimization"
run_step gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log

# --- 4. CONSTRAINED MINIMIZATION / CONJUGATE GRADIENT (CG) ---
log "Step 4: Running grompp for conjugate gradient (using gromacs-cg.mdp)"
run_step gmx grompp -v -f "$SCRIPT_DIR/gromacs-cg.mdp" -c after_em.gro -r after_em.gro -o cg -p topol.top

log "Step 4: Running mdrun for conjugate gradient"
run_step gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log

# --- 5. MOLECULAR DYNAMICS (MD) ---
log "Step 5: Running grompp for molecular dynamics (using gromacs-md.mdp)"
run_step gmx grompp -f "$SCRIPT_DIR/gromacs-md.mdp" -o md -c after_cg.gro -r after_cg.gro -p topol.top -maxwarn 1

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
run_step gmx trjconv -f system_compact.xtc -s md.tpr -o lastframe_drop.pdb -b 1000 -e 1000 <<EOF
0
EOF

# input: lastframe_drop.pdb output: filtered_renum.pdb
# Removing waters too far away from protein
log "Step 7: Removing extra waters not connected to protein"
run_step python "$SCRIPT_DIR/clean_pdb.py" "$INPUT_DIR/lastframe_drop.pdb"

log "Step 8: Post processing PDB"
run_step "$SCRIPT_DIR"/format-pdb.sh "$INPUT_DIR/filtered_renum.pdb"

log "gromacs.sh completed successfully"