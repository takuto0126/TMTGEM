#!/bin/bash

fil1=../Tohoku/em/bxyz/b14_bxyz_ts.dat
fil2=../Tohoku/em_IGRF/bxyz/b14_bxyz_ts.dat
fil3=../Tohoku/em_woa/bxyz/b14_bxyz_ts.dat
fil_ref1=Tohoku_ref/em/b14_bxyz_ts.dat
fil_ref2=Tohoku_ref/em_IGRF/b14_bxyz_ts.dat
fil_ref3=Tohoku_ref/em_woa/b14_bxyz_ts.dat

gmt begin tohoku_em pdf

# em
gmt basemap -JX20/5 -R0/5/-0.5/0.5 -Bxa1 -Bya1 -BWeSn -Y105
awk '{print($1 /60, $2)}' $fil1 | gmt plot -W1,blue
awk '{print($1 /60, $2)}' $fil_ref1 | gmt plot -Sc0.3 -W1,blue

# em_IGRF
gmt basemap -JX20/5 -R0/5/-0.5/0.5 -Bxa1 -Bya1 -BWeSn -Y-7
awk '{print($1 /60, $2)}' $fil2 | gmt plot -W1,red
awk '{print($1 /60, $2)}' $fil_ref2 | gmt plot -St0.3 -W1,red

# em_woa
gmt basemap -JX20/5 -R0/5/-0.5/0.5 -Bxa1 -Bya1 -BWeSn -Y-7
awk '{print($1 /60, $2)}' $fil_ref3 | gmt plot -Sa0.3 -W1,purple
awk '{print($1 /60, $2)}' $fil3 | gmt plot -W1,purple

gmt end show