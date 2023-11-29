# modified on 2023.11.28

grd1min=W7E28S29N44_1min.grd
xyz1min=W7E28S29N44_1min.xyz
gmt grd2xyz $grd1min -R7/28/29/44 | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

rm gmt.history

