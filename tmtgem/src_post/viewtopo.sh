# adjusted to GMT6 2021.05.21
#!/bin/bash

infile=topo.xyz
out21=topo.ps
grdf=topo.grd
WESN=-600/600/-600/600

gmt begin topo pdf
gmt xyz2grd $infile -G$grdf -R$WESN -I2
gmt makecpt -Cpolar -T-7000/7000/1

gmt basemap -JX20/20 -R${WESN} -Bxa100 -Bya100 -BWeSn
gmt grdimage $grdf  -C
gmt grdcontour $grdf -C100 -A100 -L-50/50 -W1

gmt end show
