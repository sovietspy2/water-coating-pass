# NEEDED md.tpr md.trr


gmx trjconv -f md.trr -s md.tpr -o pbc_1.xtc -pbc whole <<EOF
0
EOF

echo "Step 1 done"

gmx trjconv -f pbc_1.xtc -s md.tpr -o pbc_2.xtc -pbc cluster <<EOF
1
0
EOF

echo "Step 2 done"

gmx trjconv -f pbc_2.xtc -s md.tpr -o pbc_3.xtc -center -pbc mol -ur compact <<EOF
1
0
EOF

echo "Step 3 done"

gmx trjconv -f pbc_3.xtc -s md.tpr -o system-compact.xtc -fit progressive <<EOF
3
0
EOF

echo "Step 4 done"

gmx trjconv -f pbc_3.xtc -s md.tpr -o system_tpy.pdb -b 0 -e 0 -fit progressive <<EOF
3
0
EOF

echo "Step 5 done"

mobywat -f /home/barnabas/CLionProjects/pass2/mobywat_test/gromacs-prediction/system-compact.xtc -tpy /home/barnabas/CLionProjects/pass2/mobywat_test/gromacs-prediction/system_tpy.pdb -t A -w Auto -n 0-1000 -cls MER -m Prediction -v Diagnostic