#!/bin/bash

inp[1]=bxyz/surfcond.out
inp[2]=bxyz/bottcond.out
inp[3]=bxyz/surfcond_wocomp.out # before complementation see m_caloceancond.f90
inp[4]=bxyz/bottcond_wocomp.out # before complementation see m_caloceancond.f90
title[1]="Surf"
title[2]="Bott"
title[3]="Surf w/o complement"
title[4]="Bott w/o complement"
xx[1]=3 ; yy[1]=12
xx[2]=12 ; yy[2]=0
xx[3]=-12 ; yy[3]=-10
xx[4]=12 ; yy[4]=0
width=600
WESN=-${width}/${width}/-${width}/${width}

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
ifort tmp.f90
./a.out ## make "tmp.yzc"
########################    tmp.f end  #####################

grdf="cond.grd"
grdtopo="topo.grd"

gmt begin ocencond pdf

gmt blockmean $topofile -R$WESN -I2 | gmt surface -G$grdtopo -I2/2 -T1 -R$WESN

# surfcond
gmt makecpt -Cjet -T2.8/5.0/0.01

for i in 1 2 3 4
do
if [ -e ${inp[$i]} ] ; then
awk '{print($1,$2,$3)}' ${inp[$i]} | gmt blockmean -R$WESN -I4| gmt surface -G$grdf -I4 -T1 -R$WESN  
gmt basemap -JX8/8 -R$WESN -Bxa200+l"Distance [km]" -Bya200+l"Distance [km]" -BWeSns  -X${xx[$i]} -Y${yy[$i]} 
gmt grdimage $grdf -C  
gmt plot "$polygongmt" -W1 -Gwhite 
gmt grdcontour $grdtopo -C1000 -L-8000/-1000 -W0.2   
gmt plot "$pos5file" -L -W0.5
gmt text -JX10/10 -R0/10/0/10  -F+f12p,Helvetica+jBR -G255  <<EOF 
8   8.2  ${title[$i]} cond [S/m]
EOF
gmt colorbar -Dx8.5/0+w8/0.3 -B0.5 -C
fi
done

gmt end show

#grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN}  -X3 -Y3 > $out
rm tmp.f90 polygon2.gmt $grdtopo $grdf a.out

