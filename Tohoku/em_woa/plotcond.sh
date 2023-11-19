#!/bin/bash

inp[1]=bxyz/surfcond.out
inp[2]=bxyz/bottcond.out
inp[3]=bxyz/surfcond_wocomp.out # before complementation see m_caloceancond.f90
inp[4]=bxyz/bottcond_wocomp.out # before complementation see m_caloceancond.f90
title[1]="Surf"
title[2]="Bott"
title[3]="Surf w/o complement"
title[4]="Bott w/o complement"
out=oceancond.ps
outpdf=oceancond.pdf

width=600
WESN=-${width}/${width}/-${width}/${width}

bb=a100/a100:"distance\(km\)":WeSn
scl=15/15 ; range=0/15/0/15
scl2=17/15 ; range2=0/20/0/15

polygonfile="../mesh/polygon2.dat"
polygongmt="polygon2.gmt"
pos5file="../mesh/pos5.dat"
topofile="../mesh/topo.xyz"

cat <<EOF > tmp.f90
program tmp
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4) :: nobs
real(8),allocatable,dimension(:,:) :: xy,bxyz,z


!########################################  calculate coastline file
open(1,file="$polygonfile")
open(2,file="$polygongmt")
l=0
100 continue
read(1,*,end=99) l,npoly,iclose
read(1,*) j,cx1,cy1
write(2,'(a4)') "> -Z"
write(2,*) cy1,cx1
do i=2,npoly
read(1,*) j,cx,cy
write(2,*) cy,cx
end do
if (iclose .eq. 1)  write(2,*) cy1,cx1
goto 100
99  continue
if ( l .eq. 0 )then
write(*,*) "GEGEGE there are no lines in the input file ","$infile"
end if
close(1)
close(2)

end program tmp
EOF
gfortran tmp.f90
./a.out ## make "tmp.yzc"
########################    tmp.f end  #####################

grdf="cond.grd"
CPT="cond.cpt"
grdtopo="topo.grd"

surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN

# surfcond
makecpt -Cjet -T2.8/5.0/0.01 -V > ${CPT}

if [ -e $out ]; then
rm $out
fi

for i in 1 2 3 4
do
cat ${inp[$i]} | awk '{print($1,$2,$3)}' | surface -G$grdf -I4/4 -T1 -R$WESN
grdimage $grdf -B$bb -C${CPT} -JX${scl} -R${WESN} -V -K -X3 -Y3 >> $out
psxy "$polygongmt" -JX$scl -R$WESN  -m -K -O -V -W1,black -Gwhite >> $out
grdcontour $grdtopo -C1000 -L-8000/-1000 -W0.2,black -JX -K -O -V >> $out
psxy "$pos5file" -JX$scl -R$WESN -K -L -O -V -W0.5,black >> $out
pstext -JX$scl2 -R$range2  -W0 -G255 -K -O -V <<EOF >> $out
16.9   14.5 16 0 4 RM ${title[$i]} cond [S/m]
EOF
psscale -D15.5/6/10/0.3 -B0.1 -C$CPT -O -V >> $out
done


#grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN} -V  -X3 -Y3 > $out
rm tmp.f90 polygon2.gmt gmt.history $CPT $grdtopo $grdf a.out
ps2pdf $out $outpdf
open $outpdf &

