!------10!-------20!----
input 3d mshfile   !./em3d.msh
0homo,1cond,2model !0
homo resistivity   !100.0
output cond        !./cond_test.msh
0:elevation,1:depth!0
# of cuboid        !1
1  xminmax [km]    !  -0.2          0.2
1  yminmax [km]    !  -0.2          0.2
minmax depth [km]  !  0.8           1.1
1  rho    [Ohm.m]  ! 10.0
