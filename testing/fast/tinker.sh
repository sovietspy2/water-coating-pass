rm -rf ./tmp_t

mkdir ./tmp_t

cd ./tmp_t

wget http://files.rcsb.org/download/1R6J.pdb

cp 1R6J.pdb 1R6J_r.pdb

../../../pipeline/wdrop.sh 1R6J.pdb tinker 1R6J_r.pdb

