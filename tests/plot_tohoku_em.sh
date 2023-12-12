#!/bin/bash
# Coded on Dec 12, 2023

fil1=../Tohoku/em/bxyz/b14_bxyz_ts.dat
fil2=../Tohoku/em_IGRF/bxyz/b14_bxyz_ts.dat
fil3=../Tohoku/em_woa/bxyz/b14_bxyz_ts.dat
fil_ref1=Tohoku_ref/em/b14_bxyz_ts.dat
fil_ref2=Tohoku_ref/em_IGRF/b14_bxyz_ts.dat
fil_ref3=Tohoku_ref/em_woa/b14_bxyz_ts.dat

gmt begin tohoku_em pdf

bb="-Bxa0.5+lTime[min] -Bya0.5 -BWeSn"
# em
gmt basemap -JX20/5 -R0/5/-0.5/0.5 $bb  -Y15
awk '{print($1 /60, $2)}' $fil1     | gmt plot -Sc0.3 -W1,blue  # cal symbol
awk '{print($1 /60, $2)}' $fil_ref1 | gmt plot -W1,blue         # ref line

# em_IGRF
gmt basemap -JX20/5 -R0/5/-0.5/0.5 $bb  -Y-7
awk '{print($1 /60, $2)}' $fil2     | gmt plot -St0.3 -W1,red   # cal  symbol
awk '{print($1 /60, $2)}' $fil_ref2 | gmt plot -W1,red          # ref line

# em_woa
gmt basemap -JX20/5 -R0/5/-0.5/0.5 $bb  -Y-7
awk '{print($1 /60, $2)}' $fil3     | gmt plot  -Sa0.3 -W1,purple # cal symbol
awk '{print($1 /60, $2)}' $fil_ref3 | gmt plot -W1,purple         # ref line

gmt end show