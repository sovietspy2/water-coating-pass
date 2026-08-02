PDB="${1:-system_ref.pdb}" #defaulting to mobywat convention
TARGET="${2:-[A]}" #chain-ID target spec, e.g. [A] or [AB]

echo "Adding REMARK section to reference PDB file."
sed -i "1i\\
REMARK mobywat_reference_target ${TARGET}\\
REMARK mobywat_reference_waters Auto\\
#REMARK mobywat_min_ctol   1.0\\
#REMARK mobywat_num_ctol   4\\
#REMARK mobywat_step_ctol  0.50\\
#REMARK mobywat_min_ptol   2.50\\
#REMARK mobywat_num_ptol   1\\
#REMARK mobywat_step_ptol  0\\
#REMARK mobywat_min_stol   1.500\\
#REMARK mobywat_num_stol   1\\
#REMARK mobywat_step_stol  0
" "${PDB}"
