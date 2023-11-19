#Adjusted to 2021.05.21
#!/bin/bash
#for i in `seq -w 9000 60 9000`
#do
inputfile="./z_01_"${1}".dat"
#inputfile="./Z_LLW/z_01_00"${i}".dat"
out21="./z_01_"${1}".ps"
#out21="./Z/z_01_00"${i}".ps"
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
ifort tmp.f90
./a.out    # tmp.dat is created
#------------------------------------------------------####
WESN=10/25/30/42
grdf="eta.grd"
scl=0.9
#topogrd="topo.grd"
CPT="eta.cpt"
#makecpt -Crainbow -T-30/30/0.1 -V > ${CPT}
gmt makecpt -Cpolar -T-2/2/0.05 > ${CPT}
#xyz2grd Bathymetry.xyz -GBathymetry.grd -R0/5.488/0/3.402 -I0.014
gmt surface "tmp.dat" -G$grdf -I0.05/0.05 -T0.5 -R$WESN
#surface $etafile -G$grdeta -I2/2 -T1 -R$WESN
#surface $topofile -G$topogrd -I2/2 -T1 -R$WESN
gmt grdimage $grdf -Ba5f1WeSn -C${CPT} -Jm${scl} -K -R${WESN}  -X3 -Y1 > $out21
#grdcontour $topogrd -C200  -L-1000/0 -W1,green  -JX -V -K -O >> $out21
#psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -K -O -m -V -W2/0/0/0 >> $out21
#psxy "$pos5file" -JX$scl -R$WESN -m -K -O -V -L -W2/255/0/0 >> $out21
#psxy "obs.dat" -JX$scl -R$WESN  -G255 -K -O -P -Ss0.5c  -V -W8/0/0/0 >> $out21
#grdcontour $grdeta -C1 -A1tf6  -L0.1/5 -W0.1/0 -JX -K -V -O >> $out21
#grdcontour $grdeta -C0.5 -A0.5tf6  -L0.1/0.8 -W0.1/0 -JX -K -V -O >> $out21
#grdcontour $grdeta -C1 -A1tf6  -L-5/-0.1 -Wa0.1/0/t7_7:0 -JX -K -V -O >> $out21
gmt pscoast -R$WESN -B -Jm$scl -Dh -W0.1p/50/50/50 -G0/255/0 -K -O >> $out21
gmt psscale -Dx14.5/5/8/0.3 -B0.5  -C$CPT -K -O >> $out21
gmt psxy -R -Jm -Sc0.1 -G255/0/255 -K -O << EOF >> $out21
159.9518        41.1026
144.807595      39.0582
EOF
gmt pstext -JX17/19 -R0/17/0/19  -G255 -O <<EOF >> $out21
7 14 16 0 4 CM  Time = $minute [min]
EOF
gv $out21 &
#open $out21
rm tmp.f90 gmt.history
rm tmp.dat $grdf
rm a.out  $CPT
#done
#convert -loop 3 -delay 5 ./Z/z_01_00[5-7][0-9][0-9]0.ps tsunami.gif
