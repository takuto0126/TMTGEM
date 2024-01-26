# modified on 2023.11.28

file=`ls *.grd`
xyz1min=topo.xyz
w=`echo $file | awk -F[WESN_] '{print($2)}'`
e=`echo $file | awk -F[WESN_] '{print($3)}'`
s=`echo $file | awk -F[WESN_] '{print($4)}'`
n=`echo $file | awk -F[WESN_] '{print($5)}'`
range=$w/$e/$s/$n
echo $range
gmt grd2xyz $file -R$range | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

rm gmt.history

