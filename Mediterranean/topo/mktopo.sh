# coded by Takuto MINAMI on 2019.02.19
#!/bin/bash
# You can download topo15.grd from ftp://topex.ucsd.edu/pub/srtm15_plus/
# You can download topo30.grd from ftp://topex.ucsd.edu/pub/srtm30_plus/

W=-15
E=35
S=20
N=55

grd15sec=~/study/Backup/Topography/Scripps/topo15.grd # global topo file
grd30sec=~/study/Backup/Topography/Scripps/topo30.grd # global topo file
grd1min=~/study/Backup/Topography/Scripps/topo_19.1.grd # global topo file

range=$W/$E/$S/$N
xyz15sec=W${W}E${E}S${S}N${N}_15sec.xyz
xyz30sec=W${W}E${E}S${S}N${N}_30sec.xyz
xyz1min=W${W}E${E}S${S}N${N}_1min.xyz

grd2xyz $grd1min -R$range | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

rm gmt.history

exit


range=$W/$E/$S/$N
xyz1min=W${W}E${E}S${S}N${N}_1min.xyz
grd2xyz $grd1min -R$range  > $xyz1min



