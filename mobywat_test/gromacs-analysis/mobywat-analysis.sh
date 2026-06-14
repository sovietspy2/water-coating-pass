#!/bin/bash
# MobyWat GROMACS ANALYSIS
# Input files: md.trr, md.tpr, system_ref.pdb

gmx trjconv -f md.trr -s md.tpr -o pbc1.xtc -pbc whole << EOF
0
EOF

gmx trjconv -f pbc1.xtc -s md.tpr -o pbc2.xtc -pbc cluster << EOF
1
0
EOF

gmx trjconv -f pbc2.xtc -s md.tpr -o pbc3.xtc -center -pbc mol -ur compact << EOF
1
0
EOF

gmx confrms -one -f1 system_ref.pdb -f2 md.tpr -o fit.pdb << EOF
3
3
EOF

gmx editconf -label A -f fit.pdb -o fit.pdb

gmx trjconv -f pbc3.xtc -s fit.pdb -o system.xtc -fit progressive << EOF
3
0
EOF

gmx trjconv -f pbc3.xtc -s fit.pdb -o system_tpy.pdb -b 0 -e 0 -fit progressive << EOF
3
0
EOF

awk '
BEGIN {n=1}
$1=="ATOM" || $1=="HETATM" {
    printf "%-6s%5d%s\n", substr($0,1,6), n, substr($0,12)
    n++
    next
}
$1=="TER" {
    printf "TER\n"
    next
}
{print}
' system_ref.pdb | sponge system_ref.pdb

sed -i '1i\
REMARK mobywat_reference_target [A]\
REMARK mobywat_reference_waters Auto\
#REMARK mobywat_min_ctol   1.0\
#REMARK mobywat_num_ctol   4\
#REMARK mobywat_step_ctol  0.50\
#REMARK mobywat_min_ptol   2.50\
#REMARK mobywat_num_ptol   1\
#REMARK mobywat_step_ptol  0\
#REMARK mobywat_min_stol   1.500\
#REMARK mobywat_num_stol   1\
#REMARK mobywat_step_stol  0
' "system_ref.pdb"

echo "Done! MobyWat input files ready: system.xtc, system_tpy.pdb, fit.pdb"

mobywat -t [A] -w Auto -n 0-1000 -m Analysis -v Diagnostic