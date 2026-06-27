export TINKER_GPU=1

rm -rf ./tmp_t_gpu

mkdir ./tmp_t_gpu

cd ./tmp_t_gpu

wget http://files.rcsb.org/download/1R6J.pdb

cp 1R6J.pdb 1R6J_r.pdb

../../../pipeline/wdrop.sh 1R6J.pdb tinker SHORT 1R6J_r.pdb

