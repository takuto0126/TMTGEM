# run script for Macbook Pro 15 inch
#!/bin/sh

# intel oneapi environment
source /opt/intel/oneapi/setvars.sh

# openMP threads
export OMP_NUM_THREADS=8

src=../../tmtgem/solver
cd $src
make clean
make
cd -

rm bxyz/* exyz/* vxyz/*

time ${src}/ebfem_tsunamiEM.exe < tohoku.ctl

