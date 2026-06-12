gmx editconf -f mobywat_input.gro -o mobywat_input_tpy.pdb

gmx check -f mobywat_input.xtc

mobywat -f mobywat_input.xtc -tpy mobywat_input_tpy.pdb \
  -t A -w Auto -n 0-1000 \
  -m Prediction -cls IDa \
  -ctol 1.0 -ptol 2.5 -dmax 3.5 -mtol 1.5 -bmax 30.0 \
  -v Verbose