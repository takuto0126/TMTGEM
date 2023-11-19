# modified on 2019.02.25
# coded by Takuto MINAMI on 2018.11.13
#!/bin/bash
# You can download topo15.grd from ftp://topex.ucsd.edu/pub/srtm15_plus/
# You can download topo30.grd from ftp://topex.ucsd.edu/pub/srtm30_plus/

W=120
E=176
S=20
N=60

grd15sec=~/study/Backup/Topography/Scripps/topo15.grd # global topo file
grd30sec=~/study/Backup/Topography/Scripps/topo30.grd # global topo file
grd1min=~/study/Backup/Topography/Scripps/topo_19.1.grd # global topo file
grd2min=2min.grd

range=$W/$E/$S/$N
xyz15sec=W${W}E${E}S${S}N${N}_15sec.xyz
xyz30sec=W${W}E${E}S${S}N${N}_30sec.xyz
xyz1min=W${W}E${E}S${S}N${N}_1min.xyz
xyz2min=W${W}E${E}S${S}N${N}_2min.xyz

grd2xyz $grd1min -R$range | awk '{printf "%f %f %f \n",$1,$2,-$3}' > $xyz1min # downward positive

#grdsample $grd1min -I2m -G$grd2min -R$range -V
#grd2xyz $grd2min -R$range | awk '{printf "%f %f %f \n",$1,$2,$3}' > $xyz2min

rm gmt.history # $grd2min

