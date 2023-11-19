# Coded on 11 Nov 2018
#!/bin/bash
source /opt/intel/oneapi/setvars.sh # 2021.05.29

head="polygon"
sdir="../../tmtgem/mesh"
#ctl="mesh_500.ctl"
ctl="mesh.ctl"

# [0] coastline.exe
cd $sdir
make clean
make
cd -

./clean.sh # erase existing file

#${sdir}/msh2spherical.exe < ${ctl}
#exit

#[1]## coastline
 ${sdir}/coastline.exe < ${ctl}

#[2]## gmsh polygonki.geo
gmsh ${head}ki.geo -2 -format msh2 -bgm bgmesh.pos
gmsh ${head}ki.msh  >& /dev/null &

#[3]## mshki2ocean.f90, to extract horizontal ocean mesh
${sdir}/mshki2ocean.exe < ${ctl}
gmsh ${head}_ki.msh  >& /dev/null &

#[4]## extrude.f90, to generate ocean.msh
${sdir}/extrude.exe < ${ctl}
gmsh ocean.msh  >& /dev/null &


#[5]## mk3dgeo.f90
${sdir}/mk3dgeo.exe < ${ctl}

gmsh pre3d.geo -3 -format msh2
#gmsh pre3d.geo -3 -bgm bgmesh3d.pos

#[6]## combine3d.f90
${sdir}/combine3d.exe < ${ctl}

#[7]##
${sdir}/mkline.exe

