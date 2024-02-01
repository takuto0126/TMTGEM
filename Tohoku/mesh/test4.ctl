!------10!-------20!----
input 3d mshfile   !./em3d.msh
0homo,1cond,2model !1
input cond         !./structure/cond_G3.msh
output cond        !./cond_test.msh
0:elevation,1:depth!1
2d topo file       |../mesh_aso_A04/nakadake2dz.msh
# of cuboid        !1
1  xminmax [km]    !  -0.2          0.2
1  yminmax [km]    !  -0.2          0.2
minmax depth [km]  !  -0.3          -0.1
1  rho    [Ohm.m]  ! 10.0
