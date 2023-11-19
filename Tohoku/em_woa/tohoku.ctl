cd --------10--------20--------------------
1:err,2:fwd,3:both |2
3d em mesh file    |../mesh/em3d.msh
line infro file    |../mesh/lineinfo.dat
mesh control file  |../mesh/mesh.ctl
ocean mesh file    |../mesh/ocean.msh
start time [sec]   |0.0
time interval [sec]|5.0
end time   [sec]   |300.0
1:COMCOT,2:analytic|1
xgrid file         |../flow/layer01_x.dat
ygrid file         |../flow/layer01_y.dat
depth grid  file   |../flow/layer01.dat
header for m file  |../flow/m_01_
header for n file  |../flow/n_01_
t int COMCOT [sec] |5.0
1:vh,2:vz,3:vh+vz  |1
0:IGRF,1:FIL,2:cons|2
Fx:E,Fy:N,Fz:up[nT]|0.0            0.0            -35000.0
out folder  (bxyz) |./bxyz/
out folder  (exyz) |./exyz/
out folder  (vxyz) |./vxyz/
# of observatories |2
1:lonlath,2:xyz[km]|1
obs name           |ESA
obs position [km]  |141.35472      39.236944        0.0
obs name           |b14
obs position [km]  |144.807595     39.0582          -5.380
0:no,1:nxny,2:seafl|1
nx, ny             |200         200
time interval [sec]|60.0
1:same dpth,2:seafl|2
nx,ny-> kimsh file |../mesh/polygonki.msh
nx,ny->ki23dptr fl |../mesh/ki23dptr.dat
# of cond from top |3
air cond   [S/m]   |1.e-8
ocean cond [S/m]   |3.3
crust cond [S/m]   |0.01
1:no,2:ocecond file|2
need woa csv file  |../../woa/woa13_decav81B0_C00an01.csv

