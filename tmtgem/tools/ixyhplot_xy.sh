# Adjusted to GMT6 2021.05.21
# Coded on Nov 20, 2016
# This script make a movie composed of surface elevation, uh, and bz
#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
width=600
WESN=-${width}/${width}/-${width}/${width}
ixyhfile="./exyz/ixyh_xy2D"${1}".dat"
out21="./exyz/ixyh_xy2D"${1}".ps"
t=`echo ${1}`
echo $t
pos5file="../mesh/pos5.dat"
topofile="../mesh/topo.xyz"
polygongmt="../mesh/landpolygon.gmt"
dt=`head -8 ./easter.ctl | tail -1 | awk -F"|" '{print($2)}' `
echo "dt=" $dt "[sec]"
vhfile="./vxyz/vh"${t}".dat"
minute=` echo "scale=1; ${1} * $dt / 60 " | bc `
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
program tmp
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4) :: nobs,nx,ny,i,j
real(8),allocatable,dimension(:,:) :: xyz,xy
real(8),allocatable,dimension(:) :: ixyh,phase

open(1,file="./exyz/coord_xy2D.dat")
read(1,*) nobs
allocate(xyz(3,nobs))
do i=1,nobs
 read(1,*)xyz(1,i),xyz(2,i),xyz(3,i)
end do
close(1)
do i=1,nobs
 if (xyz(1,i) .ne. xyz(1,i+1)) then
  ny=i ; goto 101
 end if
end do
101 continue
nx = nobs /ny
write(*,*) "nx=",nx
write(*,*) "ny=",ny

allocate(ixyh(nobs),phase(nobs))
open(1,file="${ixyhfile}")
do i=1,nobs
 read(1,*) ixyh(i),phase(i)
end do
close(1)

!#[1]## amp file
open(1,file="amp.dat")
do i=1,nobs
 write(1,*) xyz(1,i),xyz(2,i),ixyh(i)
end do
close(1)

!#[2]## phase file
open(2,file="phase.dat")
do i=1,nx,3
do j=1,ny,3
 ii=(i-1)*ny + j
 if (ixyh(ii) .ge. 0.1) then
 write(2,*) xyz(1,ii),xyz(2,ii),phase(ii),0.3
 end if
end do
end do
close(2)

!#[3]## read coord for vxyz
open(1,file="./vxyz/coord_xy2D.dat")
read(1,*) nobs
allocate(xy(2,nobs))
do i=1,nobs
 read(1,*) xy(1,i),xy(2,i),z
end do
close(1)

!#[4]## read vh
open(1,file="${vhfile}")
open(2,file="vh.dat")
 do i=1,nobs
 read(1,*,end=98) a
 write(2,*) xy(1,i),xy(2,i),a
 end do
 98 continue
close(1)
close(2)
deallocate(xy)

end program tmp
EOF
########################    tmp.f end  #####################
gfortran tmp.f90
./a.out ## make "tmp.yzc"
#============ bzplot ===================
scl=15/15 ; range=0/15/0/15
bb=a400/a400:"distance\(km\)":WeSn
grdf="ixyh.grd"
grdvh="vh.grd"
CPT="bz.cpt"
gmt makecpt -Crainbow -T-7.3/15/0.01 > ${CPT}
gmt surface "amp.dat" -G$grdf -I2/2 -T1 -R$WESN
gmt surface "vh.dat" -G$grdvh -I3/3 -T0.2 -R$WESN
gmt grdimage $grdf -B$bb -C${CPT} -JX${scl} -K -R${WESN}  -X3 -Y3 > $out21
gmt grdcontour $grdvh -C10 -L10/20 -W0.2,black -JX -K -O >> $out21
gmt psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -K -O -m -W0.5,black >> $out21
gmt psxy "$pos5file" -JX$scl -R$WESN -K -L -O -W0.5/255/0/0 >> $out21
gmt psxy "phase.dat"  -R  -J -Sv0.03/0.2/0.05 -W0.1,black -G255 -K -O >> $out21
gmt psscale -D15.5/6/10/0.3 -B1 -C$CPT -G0/15 -K -O >> $out21
scl2=17/15 ; range2=0/20/0/15
gmt pstext -JX$scl2 -R$range2  -G255 -O <<EOF >> $out21
17.3   14.5 16 0 4 RM t = $minute [min]
17.3   13.8    16 0 4 RM    ${1} * ${dt} [s]
18.0 12.8 16 0 4 LM h|I_H|
18.0 12.0 16 0 4 LM [A/km]
EOF
#=========  plot end =======================
#rm ocaen.dat
rm $CPT a.out $grdf
#rm  $grdtopo
#rm vec.dat uh.grd uhxyz.dat uh.cpt eta.grd eta.cpt
rm $grdf  tmp.dat
gv $out21 &
open $out21
