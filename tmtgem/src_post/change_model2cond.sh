#!/bin/bash

# [Explanation of change_model2cond.sh ] ## by Takuto MINAMI 2022.10.25
# change_model2cond.sh can manipulate existing resistivity structures represented by tetrahedral mesh
# Input types are determined by iflag_cond(=0,1,2) while the way of manipulation is done by iflag_ele(=0,1)
#
# Input: (iflag_cond = 0,1,2)
# iflag_cond = 0 : homo geneous resistivity (the following input requires land resistivity value)
# iflag_cond = 1 : cond.msh type file (the following input resuires cond file alone)
# iflag_cond = 2 : model.dat type file (the following input requires connect file and model file)
#
# Output:
# cond.msh type (can be drawn by gmtslice_change.sh)
#
# the bottom and top of cuboid can be given by two types
# 0: elevation              : bottom and top are given by elevation
# 1: depth from the surface : depth bottom and top are given by negative values from the surface
#
# Way of manipulation (iflag_ele = 0,1)
# iflag_ele = 0 : xminmax, yminmax, and elevation minmax specify the cuboid where rho is changed
# iflag_ele = 1 : xminmax, yminmax, and minmax of depth(<0) from the surface specify the zone to be changed
#
# [Example of iflag_ele = 1]
#  The following case specify minmax depth to be -0.3 km and -0.1 km, respectively.
#  This choose area to be changed as z = surface elevation - 0.3 [km] to surface elevation - 0.1 [km]
#  See the results of Cases of [2], [4], [6]  
#============================================ part of control file of test2.ctl
#1  xminmax [km]    !  -0.2          0.2
#1  yminmax [km]    !  -0.2          0.2
#minmax depth [km]  !  -0.3          -0.1
#1  rho    [Ohm.m]  ! 10.0
#============================================

FC=gfortran

SRC=../../tmtgem/src_post
cd $SRC
make clean
make
cd -

ctlfile=test6.ctl  # choose control file from test1.ctl to test6.ctl



#[6]##  homogeneous resisticity input, depth cuboid specification
cat > test6.ctl <<EOF
!------10!-------20!----
input 3d mshfile   !./em3d.msh
0homo,1cond,2model !0
homo resistivity   !100.0
output cond        !./cond_test.msh
0:elevation,1:depth!0
# of cuboid        !1
1  xminmax [km]    !  -50.0          50.0
1  yminmax [km]    !  -50.2          50.2
minmax depth [km]  !  -0.3          -50.1
1  rho    [Ohm.m]  ! 10.0
EOF

${SRC}/change_model2cond.exe < test6.ctl
#./gmtslice_change.sh
