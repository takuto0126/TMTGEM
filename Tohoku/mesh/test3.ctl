!------10!-------20!----
input 3d mshfile   !./em3d.msh
0homo,1cond,2model !1
input cond         !./structure/cond_G3.msh
output cond        !./cond_test.msh
0:elevation,1:depth!0
# of blocks        !1
1  xminmax [km]    !  -0.20          0.20
1  yminmax [km]    !  -0.20          0.20
minmax elevation[km!   0.75          1.1
1  rho    [Ohm.m]  !  10.0
