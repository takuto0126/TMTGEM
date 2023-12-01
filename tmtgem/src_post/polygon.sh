# Coded by T. MINAMI on Jan. 10, 2014
#!/bin/bash
infile=polygon${1}.dat
intermf=polygon${1}.gmt
pos5file="pos5.dat"
topofile="topo.xyz"
out21=polygon${1}.ps
echo $infile
scale=15/15
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
    program tmp
    implicit none
    integer(4) :: l,npoly,j,i,iclose
    real(8) :: cx1,cy1,cx,cy
!###################  file read start #####################
	open(1,file="$infile")
    open(2,file="$intermf")
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
    stop
    end if
	close(1)
    close(2)
end program tmp
EOF
########################    tmp.f end  #####################
source /opt/intel/oneapi/setvars.sh
FC=ifort
$FC tmp.f90
./a.out ## make "tmp.yzc"
#############################################
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
rw=-1000
re=1000
rs=-1000
rn=1000
WESN=$rw/$re/$rs/$rn
scl=$scale
grdtopo="topo.grd"

## note that -L option enclose the polygon automatically in psxy
gmt begin polygon${1} pdf
gmt blockmean $topofile -R$WESN -I2 | gmt surface -G$grdtopo -I2 -T1 -R$WESN
gmt basemap -JX$scl -R$WESN  -Bxa100+l"distance\(km\)" -Bya100+l"distance\(km\)" -BWeSn
gmt plot $intermf  -W1 
gmt plot "$pos5file" -L -W1  # -L forces "closed"
#grdcontour $grdtopo -C1000 -A1000t -L-7500/-6500 -W0.2 -JX -O -V 
#pstext "elenum.dat"  -Jx$scl -R$WESN -K -O -V  
#pstext  -JX$scl -R$WESN -O -V -W3 <<EOF 
#2000 0 16 0 4 CM mesh
#EOF
gmt end show
rm tmp.f90
