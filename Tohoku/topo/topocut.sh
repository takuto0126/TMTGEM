#!/bin/bash


W=120
E=176
S=20
N=60
range=$W/$E/$S/$N
out=W${W}E${E}S${S}N${N}_1min_cf.grd

grd1min="/Users/minami/Dropbox/study/Backup/Topography/Scripps/topo_19.1.grd" # global topo file

gmt grdcut $grd1min -R$range -G"tmp.grd" 
gmt grdmath tmp.grd NEG = ${out}=cf
