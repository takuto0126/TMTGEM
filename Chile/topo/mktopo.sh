# modified on 2023.11.28

grd1min=W-81E-69S-53N-41_1min.grd
xyz1min=W-81E-69S-53N-41_1min.xyz
gmt grd2xyz $grd1min -R-81/-69/-53/-41 | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

rm gmt.history

