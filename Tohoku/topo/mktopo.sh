# modified on 2023.11.28

grd1min=W130E155S33N45_1min.grd
xyz1min=W130E155S33N45_1min.xyz
gmt grd2xyz $grd1min -R130/155/33/45 | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

rm gmt.history

