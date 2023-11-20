# Coded on May 18, 2016
# This script make a movie composed of surface elevation, uh, and bz
#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt

width=600
WESN=-${width}/${width}/-${width}/${width}
bxyzfile="./bxyz/bxyz_xy2D"${1}".dat"
t=`echo ${1} `
echo ${t}
vhfile="./vxyz/vh"${t}".dat"
echo $bxyzfile
out21="./bxyz/bxyz_xy2D"${1}"_"${2}".ps"
pos5file="../mesh/pos5.dat"
topofile="../mesh/topo.xyz"
polygongmt="../mesh/landpolygon.gmt"
# get dt
if [ ${2} = 1 ]; then bcomp="bx" ; fi
if [ ${2} = 2 ]; then bcomp="by" ; fi
if [ ${2} = 3 ]; then bcomp="bz" ; fi
dt=`head -8 ./tohoku.ctl | tail -1 | awk -F"|" '{print($2)}' `
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
gfortran tmp.f90
./a.out ## make "tmp.yzc"
#============ bzplot ===================
scl=15/15 ; range=0/15/0/15
bb=a100/a100:"distance\(km\)":WeSn
grdf="bz.grd"
grdvh="vh.grd"
grdtopo="topo.grd"
CPT="bz.cpt"

gmt begin "./bxyz/bxyz_xy2D"${1}"_"${2} pdf

gmt makecpt -Cpolar -T-18/18/0.05 > ${CPT}
gmt surface "vh.dat" -G$grdvh -I3/3 -T0.2 -R$WESN
gmt surface "tmp.dat" -G$grdf -I3/3 -T1 -R$WESN
gmt surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN
gmt basemap  -Bxa100+l"distance\(km\)" -Bya100+l"distance\(km\)" -BWeSn -JX${scl} -R${WESN}
gmt grdimage $grdf -C   -X3 -Y3
gmt grdcontour $grdvh -C10 -L10/20 -W0.2,black 
gmt plot "$polygongmt"  -m  -W1,black
gmt plot "$pos5file" -L -W0.5,black
gmt grdcontour $grdtopo -C1000 -L-8000/-7000 -W0.2,green

gmt colorbar -Dx15.5/6+w10/0.3 -B2 -C 
scl2=17/15 ; range2=0/20/0/15

# B14
#psxy -R -JX -Sc0.1 -G255/0/0  << EOF
#810.366308254877        622.980696018805
#EOF

gmt text -JX$scl2 -R$range2  -W0 -G255 <<EOF
16.9   14.5  t = $minute [min]
16.9   13.8  ${1} * ${dt} [s]
18.8 12.4  ${bcomp}
18.8 11.6  [nT]
EOF

gmt end show

#=========  plot end =======================

rm $CPT a.out tmp.f90 vh.dat vh.grd gmt.history tmp.dat
rm $grdtopo
rm $grdf
gv $out21 &
open $out21

