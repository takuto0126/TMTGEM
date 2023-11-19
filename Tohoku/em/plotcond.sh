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

polygongmt="../mesh/landpolygon.gmt"
pos5file="../mesh/pos5.dat"
topofile="../mesh/topo.xyz"

grdf="cond.grd"
CPT="cond.cpt"
grdtopo="topo.grd"

gmt surface $topofile -G$grdtopo -I2/2 -T1 -R$WESN

# surfcond
gmt makecpt -Cjet -T2.8/5.0/0.01 -V > ${CPT}

if [ -e $out ]; then
rm $out
fi

for i in 1 2 # 3 4
do
cat ${inp[$i]} | awk '{print($1,$2,$3)}' | gmt surface -G$grdf -I4/4 -T1 -R$WESN
gmt grdimage $grdf -B$bb -C${CPT} -JX${scl} -R${WESN} -K -X3 -Y3 >> $out
gmt psxy "$polygongmt" -JX$scl -R$WESN  -m -K -O -W1,black -Gwhite >> $out
gmt grdcontour $grdtopo -C1000 -L-8000/-1000 -W0.2,black -JX -K -O >> $out
gmt psxy "$pos5file" -JX$scl -R$WESN -K -L -O -W0.5,black >> $out
gmt pstext -JX$scl2 -R$range2  -W0 -G255 -K -O <<EOF >> $out
16.9   14.5 16 0 4 RM ${title[$i]} cond [S/m]
EOF
gmt psscale -D15.5/6/10/0.3 -B0.1 -C$CPT -O -V >> $out
done


#grdimage $grdf_y -B$bb -C${CPT} -JX${scl} -R${WESN} -V  -X3 -Y3 > $out
rm tmp.f90 gmt.history $CPT $grdtopo $grdf a.out
open $out
gv $out &

