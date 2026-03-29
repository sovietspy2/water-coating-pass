#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

command -v gmx >/dev/null || { echo "Error: gmx not found" >&2; exit 1; }
command -v python >/dev/null || { echo "Error: python not found" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_PDB="$1"
INPUT_ABS="$(cd -P "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"
INPUT_DIR="$(dirname "$INPUT_ABS")"

# --- TIMER SETUP ---
SECONDS=0
TIMEFILE="${INPUT_ABS}-mm-process-time.txt"

write_runtime() {
  local elapsed="$SECONDS"
  # Write runtime (in seconds) to the requested file
  printf 'runtime for mm is %s seconds\n' "$elapsed" > "$TIMEFILE"
}
trap write_runtime EXIT

cd "$INPUT_DIR"

# --- 1. SETUP ---
# Generate topology (topol.top) and initial structure (conf.gro)
gmx pdb2gmx -water tip3p -ff amber99sb-ildn -ignh -f "$INPUT_ABS" -o conf.gro

# --- 2. PREPARATION FOR PSEUDO-PBC ---
# Center system and place it in a large cubic box (1000 nm) to satisfy the 333.3 nm cutoffs
# This replaces the original gmx editconf/solvate steps.
gmx editconf -f conf.gro -o conf_largebox.gro -c -d 0 -bt cubic -box 1000

# --- 3. ENERGY MINIMIZATION (EM) ---
# Use the pseudo-PBC settings (large cutoffs) in gromacs-st.mdp
gmx grompp -v -f "$SCRIPT_DIR/gromacs-st.mdp" -c conf_largebox.gro -r conf_largebox.gro -o em -p topol.top -maxwarn 2
gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log

# --- 4. CONSTRAINED MINIMIZATION / CONJUGATE GRADIENT (CG) ---
# Use the structure from the EM, still in the large box.
# **Must use pseudo-PBC settings in gromacs-cg.mdp**
gmx grompp -v -f "$SCRIPT_DIR/gromacs-cg.mdp" -c after_em.gro -r after_em.gro -o cg -p topol.top
gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log

# --- 5. MOLECULAR DYNAMICS (MD) ---
# Use the structure from the CG, still in the large box.
# **Must use pseudo-PBC settings in gromacs-md.mdp**
gmx grompp -f "$SCRIPT_DIR/gromacs-md.mdp" -o md -c after_cg.gro -r after_cg.gro -p topol.top -maxwarn 1
gmx mdrun -v -s md -o md.trr -c after_md.gro -g md.log

# --- 6. POST-PROCESSING (Cleaning up the trajectory) ---
# The large box is unnecessary for analysis.
# We correct for the artificial PBC due to the 1000 nm box.

# 6a. Center the whole peptide+water molecule to remove jumps caused by the large box
gmx trjconv -f md.trr -s md.tpr -o pbc_whole.xtc -pbc whole <<EOF
0
EOF
# 6b. Center the system and remove the periodic box for a compact view.
# '1' is the Solvent group (water), '0' is System.
gmx trjconv -f pbc_whole.xtc -s md.tpr -o system_compact.xtc -center -pbc mol -ur compact -boxcenter zero <<EOF
1
0
EOF

# Optional: Create a final frame pdb
gmx trjconv -f system_compact.xtc -s md.tpr -o lastframe_drop.pdb -b 1000 -e 1000 <<EOF
0
EOF

#removing extra waters that are not connected to the protein
python "$SCRIPT_DIR/clean_pdb.py" "$INPUT_DIR/lastframe_drop.pdb"

echo "gromacs done"