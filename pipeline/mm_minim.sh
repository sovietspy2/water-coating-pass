#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.pdb>" >&2
  exit 1
fi

INPUT_PDB="$1"

# --- 1. SETUP ---
# Generate topology (topol.top) and initial structure (conf.gro)
gmx pdb2gmx -water tip3p -ff amber99sb-ildn -ignh -f  "$INPUT_PDB" -o conf.gro

# --- 2. PREPARATION FOR PSEUDO-PBC ---
# Center system and place it in a large cubic box (1000 nm) to satisfy the 333.3 nm cutoffs
# This replaces the original gmx editconf/solvate steps.
gmx editconf -f conf.gro -o conf_largebox.gro -c -d 0 -bt cubic -box 1000

# --- 3. ENERGY MINIMIZATION (EM) ---
# Use the pseudo-PBC settings (large cutoffs) in st.mdp
gmx grompp -v -f st.mdp -c conf_largebox.gro -r conf_largebox.gro -o em -p topol.top -maxwarn 2
gmx mdrun -v -s em -o em.trr -c after_em.gro -g em.log

# --- 4. CONSTRAINED MINIMIZATION / CONJUGATE GRADIENT (CG) ---
# Use the structure from the EM, still in the large box.
# **Must use pseudo-PBC settings in cg.mdp**
gmx grompp -v -f cg.mdp -c after_em.gro -r after_em.gro -o cg -p topol.top
gmx mdrun -v -s cg -o cg.trr -c after_cg.gro -g cg.log

# --- 5. MOLECULAR DYNAMICS (MD) ---
# Use the structure from the CG, still in the large box.
# **Must use pseudo-PBC settings in md.mdp**
gmx grompp -f md.mdp -o md -c after_cg.gro -r after_cg.gro -p topol.top -maxwarn 1
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

python clean_pdb.py