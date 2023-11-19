#!/bin/bash

inp=bxyz/geomag.out
out=geomag.ps
outpdf=geomag.pdf

width=1200
WESN=-${width}/${width}/-${width}/${width}
#WESN=-1200/0/-400/800

bb=a400/a400:"distance\(km\)":WeSn
scl=15/15 ; range=0/15/0/15
scl2=17/15 ; range2=0/20/0/15

polygongmt="../mesh_large/landpolygon.gmt"
pos5file="../mesh_large/pos5.dat"

grdf_x="igrf_x.grd"
grdf_y="igrf_y.grd"
grdf_z="igrf_z.grd"

CPT="igrf.cpt"

cat $inp | awk '{print($1,$2,$3)}' | surface -G$grdf_x -I3/3 -T0.2 -R$WESN
cat $inp | awk '{print($1,$2,$4)}' | surface -G$grdf_y -I3/3 -T0.2 -R$WESN
cat $inp | awk '{print($1,$2,$5)}' | surface -G$grdf_z -I3/3 -T0.2 -R$WESN

# Fx
makecpt -Crainbow -T-5000/-2000/100 -V > ${CPT}
grdimage $grdf_x -B$bb -C${CPT} -JX${scl} -R${WESN} -V -K -X3 -Y3 > $out
psxy "$polygongmt" -JX$scl -R$WESN  -m -K -O -V -W1,black -Gwhite >> $out
psxy "$pos5file" -JX$scl -R$WESN -K -L -O -V -W0.5,black >> $out
pstext -JX$scl2 -R$range2  -W0 -G255 -K -O -V <<EOF >> $out
16.9   14.5 16 0 4 RM Fx [nT]
EOF
psscale -D15.5/6/10/0.3 -B1000 -C$CPT -O -V >> $out

# Fy
makecpt -Crainbow -T24000/34000/500 -V > ${CPT}
grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN} -V -K -X3 -Y3 >> $out
psxy "$polygongmt" -JX$scl -R$WESN  -m -K -O -V -W1,black -Gwhite >> $out
psxy "$pos5file" -JX$scl -R$WESN -K -L -O -V -W0.5,black >> $out
pstext -JX$scl2 -R$range2  -W0 -G255 -K -O -V <<EOF >> $out
16.9   14.5 16 0 4 RM Fy [nT]
EOF
psscale -D15.5/6/10/0.3 -B1000 -C$CPT -O -V >> $out

# Fz
makecpt -Crainbow -T-45000/-30000/500 -V > ${CPT}
grdimage $grdf_z -B$bb -C${CPT} -JX${scl} -R${WESN} -V -K -X3 -Y3 >> $out
psxy "$polygongmt" -JX$scl -R$WESN  -m -K -O -V -W1,black -Gwhite >> $out
psxy "$pos5file" -JX$scl -R$WESN -K -L -O -V -W0.5,black >> $out
pstext -JX$scl2 -R$range2  -W0 -G255 -K -O -V <<EOF >> $out
16.9   14.5 16 0 4 RM Fz [nT]
EOF
psscale -D15.5/6/10/0.3 -B1000 -C$CPT -O -V >> $out


#grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN} -V  -X3 -Y3 > $out
rm tmp.f90 gmt.history $grdf_x $grdf_y $grdf_z a.out $CPT
ps2pdf $out $outpdf
open $outpdf &

