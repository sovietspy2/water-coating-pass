#!/bin/bash
# MobyWat GROMACS Trajectory Preprocessor
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

gmx confrms -one -f1 system_reference.pdb -f2 md.tpr -o fit.pdb << EOF
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

echo "Done! MobyWat input files ready: system.xtc, system_tpy.pdb, fit.pdb"

mobywat -tpy system_tpy.pdb -t [A] -w Auto -n 0-100 -m Analysis -v Diagnostic