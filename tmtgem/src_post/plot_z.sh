#!/bin/bash
#for i in `seq -w 9000 60 9000`
#do
inputfile="./z_01_"${1}".dat"
#inputfile="./Z_LLW/z_01_00"${i}".dat"
nx=`wc layer01_x.dat | awk -F" "  '{printf($1)}' `
ny=`wc layer01_y.dat | awk -F" "  '{printf($1)}' `
echo $inputfile
echo "nx=" $nx "ny=" $ny
minute=`echo "scale=2; "${1}"/60" | bc | sed -e 's/^\./0./g' `    # when dt = 1sec
#echo $minute
#------------------------------------------------------
cat <<EOF > tmp.f90
program tmp
implicit none
integer(4),parameter :: nx=${nx},ny=${ny}
real(8)   :: eta(nx,ny), x(nx), y(ny)
character(50) :: nxfile="layer01_x.dat"
character(50) :: nyfile="layer01_y.dat"
character(50) :: zfile="$inputfile"
character(50) :: outfile="tmp.dat"
integer(4) :: i,j
!#[1]##  Read coordinates
open(1,file=nxfile)
do i=1,nx
  read(1,*) x(i)
end do
close(1)
open(2,file=nyfile)
do i=1,ny
  read(2,*)y(i)
end do
close(2)
!#[2]##  Read sea surface displacement
!#
open(3,file=zfile)
do i=1,ny
  read(3,'(15f9.4)') (eta(j,i),j=1,nx)
end do
close(3)
!#[3]## Output the z data with coordinates
open(4,file=outfile)
write(4, '(3g15.7)' ) ((x(i),y(j),eta(i,j), i=1,nx),j=1,ny)
close(4)
write(*,*) "tmp.f end!!"
end program tmp
EOF
#------------------------------------------------------####
source /opt/intel/oneapi/setvars.sh
FC=ifort
$FC tmp.f90
./a.out    # tmp.dat is created
#------------------------------------------------------####
w=`head -75 comcot.ctl  |tail -1| awk -F":" '{print($2)}'| sed 's/ *$//g'| sed 's/^ *//g'`
e=`head -76 comcot.ctl  |tail -1| awk -F":" '{print($2)}'| sed 's/ *$//g'| sed 's/^ *//g'`
s=`head -77 comcot.ctl  |tail -1| awk -F":" '{print($2)}'| sed 's/ *$//g'| sed 's/^ *//g'`
n=`head -78 comcot.ctl  |tail -1| awk -F":" '{print($2)}'| sed 's/ *$//g'| sed 's/^ *//g'`
WESN=$w/$e/$s/$n
echo $WESN

grdf="eta.grd"

gmt begin z_01_${1} pdf

gmt makecpt -Cpolar -T-2/2/0.05
#xyz2grd Bathymetry.xyz -GBathymetry.grd -R0/5.488/0/3.402 -I0.014
gmt blockmean "tmp.dat" -I0.02 -R$WESN | gmt surface -G$grdf -I0.02 -T0.5 -R$WESN
#surface $etafile -G$grdeta -I2/2 -T1 -R$WESN
#surface $topofile -G$topogrd -I2/2 -T1 -R$WESN
gmt basemap -Bxa5+l"Longitude" -Bya5+l"Latitude" -BWeSn -R${WESN} -JM13
gmt grdimage $grdf  -C 
#grdcontour $topogrd -C1000  -L-7100/-6900 -W0.1/0/255/0  -JX 
#psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -m -W2/0/0/0 
#psxy "$pos5file" -JX$scl -R$WESN -m -L -W2/255/0/0 
#psxy "obs.dat" -JX$scl -R$WESN  -G255 -P -Ss0.5c  -W8/0/0/0 
gmt coast -Dh -W0.1 -Ggreen 
gmt colorbar -Dx13.5/0+w8/0.3 -B0.5  -C 
gmt plot -G255/0/255 << EOF 
159.9518        41.1026
144.807595      39.0582
EOF
gmt text -JX17/19 -R0/17/0/19  -F+f16p,Helvetica+jRB -G255 <<EOF 
13 0.2 Time = $minute [min]
EOF

gmt end show

rm tmp.f90 
rm tmp.dat $grdf
rm a.out  
