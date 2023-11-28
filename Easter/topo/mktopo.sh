# modified on 2023.11.28

grd1min=W-121E-97S-38N-16_1min.grd
xyz1min=W-121E-97S-38N-16_1min.xyz
gmt grd2xyz $grd1min -R-121/-97/-38/-16 | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

rm gmt.history

