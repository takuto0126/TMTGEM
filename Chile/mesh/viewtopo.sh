# Adjusted to GMT6 2021.05.24
#!/bin/bash

infile=topo.xyz
out21=topo.ps
grdf=topo.grd
WESN=-2000/2000/-3000/3000
CPT=topo.cpt

gmt xyz2grd $infile -G$grdf -R$WESN -I2
gmt makecpt -Cpolar -T-7000/7000/1 > ${CPT}

gmt grdimage $grdf -Ba1000/a1000WeSn -C${CPT} -JX20/20 -R${WESN}  -X3 -Y1 > $out21

#open $out21 &
gv $out21 &
