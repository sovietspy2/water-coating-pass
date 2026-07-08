rm -rf ./tmp_g

mkdir ./tmp_g

cd ./tmp_g

wget http://files.rcsb.org/download/1R6J.pdb

cp 1R6J.pdb 1R6J_r.pdb

../../../pipeline/wdrop.sh 1R6J.pdb gromacs --layers 5 --intermediate-md-ns 0 --final-md-ns 0.002 1R6J_r.pdb

