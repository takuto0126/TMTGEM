# Coded on May 18, 2016
# This script make a movie composed of surface elevation, uh, and bz
#!/bin/bash
export PATH=$PATH:/usr/lib/gmt/bin
export GMTHOME=/usr/lib/gmt
width=400
WESN=-${width}/${width}/-${width}/${width}
# get dt

pos5file="../mesh/pos5.dat"
topofile="../mesh/topo.xyz"
polygonfile="../mesh/polygon2.dat"
polygongmt="polygon2.gmt"
dt=5  # time step of TMTGEM simulation
echo "dt=" $dt "[sec]"
out21="./bh_xy2D.ps"
tstep=$1

scl=15/15 ; range=0/15/0/15
scl2=17/15 ; range2=0/20/0/15
scl3=20/19 ; range3=0/20/0/15
gmt pstext -JX$scl2 -R$range2  -W0 -G255 -K <<EOF > $out21
16.9   14.5 16 0 4 RM  [min]
EOF

echo "tstep=" ${tstep}
bxyzfile="./bxyz/bxyz_xy2D"${tstep}".dat"
t=`echo ${tstep} | cut -c 3-6`
echo ${t}
vhfile="./vxyz/vh"${t}".dat"
zfile="./eta/z"${i}".dat"
echo $bxyzfile
minute=` echo "scale=1; ${tstep} * $dt / 60 " | bc `
#######################  tmp2.f start  ###############
cat <<EOF > tmp.f90
program tmp
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4) :: nobs,nx,ny
real(8),allocatable,dimension(:,:) :: xy,bxyz,z
real(8) :: b,phase,r2d
r2d=180./(4.*datan(1.d0))
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

!# amp file
open(1,file="tmp.dat")
do i=1,nobs
 b = sqrt(bxyz(1,i)**2.d0 + bxyz(2,i)**2.d0)
 if ( b .ge. 5.d0) b=4.999
 write(1,*) xy(1,i),xy(2,i),b
end do
close(1)


!#[2]## phase file
open(2,file="phase.dat")
nx=0
do i=1,nobs
 nx = nx + 1
 if (xy(1,i) .ne. xy(1,i+1)) goto 110
end do
110 continue
ny = nobs/nx
write(*,*) "nx,ny=",nx,ny

isamp=4
!#if ( ${i} .le. 1 ) isamp = 4
do i=1,nx,isamp
do j=1,ny,isamp
  ii=(i-1)*ny + j
 b = sqrt(bxyz(1,ii)**2.d0 + bxyz(2,ii)**2.d0)
 if ( b .gt. 0.15 ) then
  phase = atan2(bxyz(2,ii),bxyz(1,ii))*r2d
  write(2,'(4g15.7)') xy(1,ii),xy(2,ii),phase,0.2
 end if
 210 continue
end do
end do
close(2)

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
########################    tmp.f end  #####################
gfortran tmp.f90
./a.out ## make "tmp.yzc"
#============ bzplot ===================
#bb=a400/a400:"distance\(km\)":WeSn
bb=a200/a200WeSn
ANOT_FONT_SIZE=10
grdf="bz.grd"
grdz="z.grd"
grdtopo="topo.grd"
CPT="bz.cpt"
gmt gmtset FONT_ANNOT_PRIMARY 12p
gmt gmtset MAP_FRAME_PEN 1p
gmt makecpt -Crainbow -T-1.6/6/0.05 > ${CPT}
#makecpt -Cgray -T-1.0/1.0/0.05 -V > ${CPT}
#gmt surface "$zfile" -G$grdz -I2/2 -T0.2 -R$WESN
gmt blockmean "tmp.dat" -R${WESN} -I2/2 > "tmp2.dat"
gmt surface "tmp.dat" -G$grdf -I2/2 -T1 -R$WESN
#surface "vh.dat" -G$grdvh -I10/10 -T0.2 -R$WESN
#surface "tmp.dat" -G$grdf -I10/10 -T1 -R$WESN
#gmt surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN
gmt grdimage $grdf -B$bb -C${CPT} -JX${scl} -K -O -R${WESN}  >> $out21
gmt psxy "phase.dat"  -R  -J -Sv0.02/0.1/0.07 -W0.2 -K -O >> $out21
#gmt grdcontour $grdz -C0.2 -L0.2/0.4 -W0.2/0/255/0 -JX -K -O >> $out21
gmt psxy "$polygongmt" -JX$scl -R$WESN  -B$ant -K -O -m -W0.5,black >> $out21 ## topography
#gmt psxy "$pos5file" -JX$scl -R$WESN -K -L -O -W0.5/255/0/0 >> $out21
#gmt grdcontour $grdtopo -C1000 -L-8000/-7000 -W0.2/0/255/0 -JX -K -O >> $out21
gmt psxy -R -JX -Ss0.2 -G200/0/0 -K -O << EOF >> $out21
810.366308254877        622.980696018805
EOF
gmt pstext -JX$scl2 -R$range2  -W0 -G255 -K -O <<EOF >> $out21
16.9   14.5 14 0 4 RM $minute [min]
EOF

gmt pstext -JX$scl3 -R$range3  -G255 -K -O <<EOF >> $out21
16.2 10.2 25 0 4 CM bh
16.2 9.4 16 0 4 CM [nT]
EOF
gmt psscale -Dx26/6/10/0.3 -B1 -G0/6 -C$CPT -O -X-10 >> $out21

#=========  plot end =======================
#rm ocaen.dat
rm $CPT a.out
rm $polygongmt #$grdtopo
rm $grdf  #uhxyz.dat
convert -density 300 "./bh_xy2D.ps" out.pdf
#gv $out21 &
open out.pdf
