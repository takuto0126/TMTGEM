!------10!-------20!----
input 3d mshfile   !./em3d.msh
0homo,1cond,2model !0
homo resistivity   !100.0
output cond        !./cond_test.msh
0:elevation,1:depth!0
# of cuboid        !1
1  xminmax [km]    !  -50.0          50.0
1  yminmax [km]    !  -50.2          50.2
minmax depth [km]  !  -50.3           0.1
1  rho    [Ohm.m]  ! 10.0
