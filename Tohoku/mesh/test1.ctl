!------10!-------20!----
input 3d mshfile   !./em3d.msh
0homo,1cond,2model !2
model connect file !./result_inv/model_connect.dat
input model        !./result_inv/model03.dat
output cond        !./cond_test.msh
0:elevation,1:depth!0
# of cuboid        !1
1  xminmax [km]    !  -0.2          0.2
1  yminmax [km]    !  -0.2          0.2
minmax elevation[km!  0.75          1.1
1  rho    [Ohm.m]  ! 10.0
