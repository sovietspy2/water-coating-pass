# Needed system_ref.pdb md.trr and md.tpr


gmx confrms -one -f1 system_ref.pdb -f2 md.tpr -o fit.pdb <<EOF
3
3
EOF

gmx editconf -label A -f fit.pdb -o fit.pdb

gmx trjconv -f md.trr -s md.tpr -o pbc_1.xtc -pbc whole <<EOF
0
EOF

gmx trjconv -f pbc_1.xtc -s md.tpr -o pbc_2.xtc -pbc cluster <<EOF
1
0
EOF

gmx trjconv -f pbc_2.xtc -s md.tpr -o pbc_3.xtc -center -pbc mol -ur compact <<EOF
1
0
EOF

gmx trjconv -f pbc_3.xtc -s fit.pdb -o system.xtc -fit progressive <<EOF
3
0
EOF

gmx trjconv -f pbc_3.xtc -s fit.pdb -o system_tpy.pdb -b 0 -e 0 -fit progressive <<EOF
3
0
EOF


mobywat -tpy system_tpy.pdb -t [A] -w Auto -n 0-100 -m Analysis -v Diagnostic