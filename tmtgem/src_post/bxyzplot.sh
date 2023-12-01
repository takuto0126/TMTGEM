# Coded on May 18, 2016
# This script make a movie composed of surface elevation, uh, and bz
#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt

bxyzfile="./bxyz/bxyz_xy2D"${1}".dat"
t=`echo ${1} `
echo ${t}
vhfile="./vxyz/vh"${t}".dat"
echo $bxyzfile
pos5file="../mesh/pos5.dat"
topofile="../mesh/topo.xyz"
polygongmt="../mesh/landpolygon.gmt"
# get dt
if [ ${2} = 1 ]; then bcomp="Bx" ; fi
if [ ${2} = 2 ]; then bcomp="By" ; fi
if [ ${2} = 3 ]; then bcomp="Bz" ; fi
dt=`head -8 ./tmtgem.ctl | tail -1 | awk -F"|" '{print($2)}' `
echo "dt=" $dt "[sec]"
minute=` echo "scale=1; ${1} * $dt / 60 " | bc `
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
program tmp
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4) :: nobs
real(8),allocatable,dimension(:,:) :: xy,bxyz,z

!#[1]## read coord for vxyz
open(1,file="./vxyz/coord_xy2D.dat")
read(1,*) nobs
allocate(xy(2,nobs))
do i=1,nobs
 read(1,*) xy(1,i),xy(2,i),z
end do
close(1)

!#[2]## read vh
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

!#[3]##
open(1,file="./bxyz/coord_xy2D.dat")
read(1,*) nobs
allocate(xy(2,nobs))
do i=1,nobs
 read(1,*)xy(1,i),xy(2,i),z
end do
close(1)

allocate(bxyz(3,nobs))
open(1,file="${bxyzfile}")
do i=1,nobs
 read(1,*) bxyz(1,i),bxyz(2,i),bxyz(3,i)
end do
close(1)

open(1,file="tmp.dat")
do i=1,nobs
 write(1,*) xy(1,i),xy(2,i),bxyz(${2},i)
end do
close(1)

end program tmp
EOF
########################    tmp.f end  #####################
source /opt/intel/oneapi/setvars.sh 
FC=ifort
$FC tmp.f90
./a.out ## make "tmp.yzc"

#============ bzplot ===================

grdf="bz.grd"
grdvh="vh.grd"
grdtopo="topo.grd"

width=600
WESN=-${width}/${width}/-${width}/${width}

gmt begin "./bxyz/bxyz_xy2D"${1}"_"${2} pdf

gmt gmtset FONT_ANNOT_PRIMARY="14p,Helvetica,black"
gmt gmtset FONT_LABEL="18p,Helvetica,black"

gmt makecpt -Cpolar -T-18/18/0.05
gmt blockmean "vh.dat"  -R$WESN -I3 |gmt surface -G$grdvh   -I3 -T0.2  -R$WESN
gmt blockmean "tmp.dat" -R$WESN -I3 |gmt surface -G$grdf    -I3 -T1    -R$WESN
gmt blockmean $topofile -R$WESN -I3 |gmt surface -G$grdtopo -I2 -T1    -R$WESN
gmt basemap  -Bxa100+l"Distance [km]" -Bya100+l"Distance [km]" -BWeSn -JX15/15 -R${WESN} -X3 -Y3
gmt grdimage $grdf -C   
gmt grdcontour $grdvh -C10 -L10/20 -W0.2
gmt plot "$polygongmt"  -W1 -L              # coastline
gmt plot "$pos5file" -L -W0.5
gmt grdcontour $grdtopo -C1000 -L-8000/-7000 -W0.2,green # bathymetry contour

gmt colorbar -Dx15.5/0+w10/0.3 -B2 -C 

gmt text -JX17/15 -R0/17/0/15 -F+f16p,Helvetica+jRM -Gwhite <<EOF
14.8  14.6  t = $minute [min]
14.8  13.9  ${dt} * ${1} [s]
16.4 11.6  ${bcomp}
16.5 10.8  [nT]
EOF

gmt end show

#=========  plot end =======================

rm $CPT a.out tmp.f90 vh.dat vh.grd tmp.dat
rm $grdtopo
rm $grdf

