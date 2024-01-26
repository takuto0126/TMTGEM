#!/bin/bash

fil1=~/jenkins/TMTGEM/Easter/em/bxyz/IPM_bxyz_ts.dat
fil2=Easter_ref/IPM_bxyz_ts.dat

gmt begin easter_em pdf

gmt basemap -JX20/10 -R0/5/-1/1 -Bxa1 -Bya1 -BWeSn
awk '{print($1 /60, $2)}' $fil1 | gmt plot -W1,blue
awk '{print($1 /60, $2)}' $fil2 | gmt plot -W1,red

gmt end show