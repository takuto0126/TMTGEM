#Adjusted to GMT6 2021.05.24
#!/bin/bash

inp=bxyz/geomag.out
out=geomag.ps
outpdf=geomag.pdf

width=800
WESN=-${width}/${width}/-${width}/${width}
#WESN=-1200/0/-400/800

bb=a400/a400:"distance\(km\)":WeSn
scl=15/15 ; range=0/15/0/15
scl2=17/15 ; range2=0/20/0/15

polygongmt="../mesh/landpolygon.gmt"
pos5file="../mesh/pos5.dat"

grdf_x="igrf_x.grd"
grdf_y="igrf_y.grd"
grdf_z="igrf_z.grd"

CPT="igrf.cpt"

cat $inp | awk '{print($1,$2,$3)}' | gmt surface -G$grdf_x -I10/10 -T1 -R$WESN
cat $inp | awk '{print($1,$2,$4)}' | gmt surface -G$grdf_y -I10/10 -T1 -R$WESN
cat $inp | awk '{print($1,$2,$5)}' | gmt surface -G$grdf_z -I10/10 -T1 -R$WESN
Fxmax=`cat $inp | awk '{if(m<$3) m=$3} END{print m + 1000 }'`
Fxmin=`cat $inp | awk 'BEGIN{m=100000}{if(m>$3) m=$3} END{print m - 1000}'`
Fymax=`cat $inp | awk '{if(m<$4) m=$4} END{print m + 1000 }'`
Fymin=`cat $inp | awk 'BEGIN{m=100000}{if(m>$4) m=$4} END{print m - 1000}'`
Fzmax=`cat $inp | awk '{if(m<$5) m=$5} END{print m + 1000 }'`
Fzmin=`cat $inp | awk 'BEGIN{m=100000}{if(m>$5) m=$5} END{print m - 1000}'`

# Fx
gmt makecpt -Crainbow -T${Fxmin}/${Fxmax}/100 -V > ${CPT}
gmt grdimage $grdf_x -B$bb -C${CPT} -JX${scl} -R${WESN}  -K -X3 -Y3 > $out
gmt psxy "$polygongmt" -JX$scl -R$WESN  -K -O -W1,black -m -Gwhite >> $out
gmt psxy "$pos5file" -JX$scl -R$WESN -K -L -O -W0.5,black >> $out
gmt pstext -JX$scl2 -R$range2  -W0 -G255 -K -O <<EOF >> $out
16.9   14.5 16 0 4 RM Fx [nT]
EOF
gmt psscale -Dx15.5/6/10/0.3 -B2000 -C$CPT -O -V >> $out

# Fy
gmt makecpt -Crainbow -T${Fymin}/${Fymax}/1000 > ${CPT}
gmt grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN} -K -X3 -Y3 >> $out
gmt psxy "$polygongmt" -JX$scl -R$WESN -m -K -O -W1,black -Gwhite >> $out
gmt psxy "$pos5file" -JX$scl -R$WESN -K -L -O -W0.5,black >> $out
gmt pstext -JX$scl2 -R$range2  -W0 -G255 -K -O <<EOF >> $out
16.9   14.5 16 0 4 RM Fy [nT]
EOF
gmt psscale -Dx15.5/6/10/0.3 -B2000 -C$CPT -O >> $out

# Fz
gmt makecpt -Crainbow -T${Fzmin}/${Fzmax}/1000 > ${CPT}
gmt grdimage $grdf_z -B$bb -C${CPT} -JX${scl} -R${WESN} -K -X3 -Y3 >> $out
gmt psxy "$polygongmt" -JX$scl -R$WESN -m -K -O -W1,black -Gwhite >> $out
gmt psxy "$pos5file" -JX$scl -R$WESN -K -L -O -W0.5,black >> $out
gmt pstext -JX$scl2 -R$range2  -W0 -G255 -K -O <<EOF >> $out
16.9   14.5 16 0 4 RM Fz [nT]
EOF
gmt psscale -Dx15.5/6/10/0.3 -B2000 -C$CPT -O >> $out


#grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN} -V  -X3 -Y3 > $out
rm  gmt.history $grdf_x $grdf_y $grdf_z $CPT
#ps2pdf $out $outpdf
gv $out &
#open $out

