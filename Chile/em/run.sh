# run script for Macbook Pro 15 inch
#!/bin/sh
source /opt/intel/bin/compilervars.sh intel64
source /opt/intel/bin/debuggervars.sh intel64

OMP_NUM_THREADS=4


src=../../tmtgem/solver
cd $src
#make clean
make
cd -

rm bxyz/* exyz/* vxyz/*

time ${src}/ebfem_tsunamiEM.exe < chile.ctl

