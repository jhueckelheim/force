#!/bin/sh
module use /soft/modulefiles
module load gcc/9.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hueckelh/mpfr-4.1.0/install/lib

#MPFR_HOME=/home/hueckelh/mpfr-4.1.0/install make --always-make

export OMP_PLACES=sockets
export OMP_NUM_THREADS=28;

for i in 1024 2048 4096 8192 16384 32768 65536 131072 262144
do
  ./exec_test_mxm $i
done

