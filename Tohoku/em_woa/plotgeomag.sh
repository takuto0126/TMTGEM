#!/bin/bash

inp=bxyz/geomag.out

width=600
WESN=-${width}/${width}/-${width}/${width}


polygonfile="../mesh/polygon2.dat"
polygongmt="polygon2.gmt"
pos5file="../mesh/pos5.dat"

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

grdf[1]="igrf_x.grd"
grdf[2]="igrf_y.grd"
grdf[3]="igrf_z.grd"

CPT="igrf.cpt"

gmt begin geomag pdf

awk '{print($1,$2,$3)}' $inp | gmt surface -G${grdf[1]} -I3/3 -T0.2 -R$WESN
awk '{print($1,$2,$4)}' $inp | gmt surface -G${grdf[2]} -I3/3 -T0.2 -R$WESN
awk '{print($1,$2,$5)}' $inp | gmt surface -G${grdf[3]} -I3/3 -T0.2 -R$WESN

xx[1]=3   ; Comp[1]=Fx ; TR[1]=-20000/20000/100
xx[2]=14  ; Comp[2]=Fy ; TR[2]=0/50000/100
xx[3]=14  ; Comp[3]=Fz ; TR[3]=-40000/0/100

for i in 1 2 3
do

# Fx
gmt makecpt -Crainbow -T${TR[$i]}
gmt basemap -Bxa100 -Bya100 -BWeSn -JX10/10 -R$WESN -X${xx[$i]}
gmt grdimage ${grdf[$i]} -C   
gmt plot "$polygongmt" -W1 -Gwhite 
gmt plot "$pos5file" -L -W0.5
gmt text -JX$17/15 -R0/17/0/15 -G255  <<EOF 
10.3   10.5  ${Comp[$i]} [nT]
EOF
gmt colorbar -Dx10.2/0+w10/0.3 -B10000 -C 

done

gmt end show

rm tmp.f90 $polygongmt gmt.history $grdf_x $grdf_y $grdf_z a.out $CPT

