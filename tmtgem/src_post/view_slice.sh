# coded on 2017.12.20
# to draw a slice
#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

igendat=1 # 0: no dat file generation, 1: dat file generation
bl=-1.0
br=1.0
bb=0.25
bt=1.5
SL1="WE"
A=0.0
B=1.0
D=0.0
SL2="SN"
A2=1.0
B2=0.0
D2=0.0

ID=test
gmt gmtset FONT_ANNOT_PRIMARY="7p,Helvetica.black"

#[1]## generate control file
# slice plane Ax + By + Cz + D = 0
cat << EOF > tmp.ctl
##--------------- control file for slice.exe starts here!! ----------------------
## 2020.12.10 Takuto Minami
##
## Note that lines starting with "##" work for comments ---
##
##             20->|
mshfile            !../mesh_aso_A04/nakadake3d.msh
## You should chose either of two types of resistivity structure file format,
## cond.msh type or model**.dat type
## 0 for cond.msh type, 1 for model**.dat type
##  case for cond.msh type
0:cond,1:model     !0
ncondfile          !1
polygonhead        !$ID
condfile           !./cond_test.msh
## case for model**.dat type
##0:cond,1:model     !1
##model connect      !./result_inv/model_connect.dat
##nmodelfile         !1
## polygonhead: header of output file
##polygonhead        !$ID
##modelfile          !./result_inv/model03.dat
## ibound 0: no bound, 1:lower, 2: upper bound
## If you want not to show some resistivity values larger / smaller than some threshold
## you can set by ibound =1 for lower bound, and 2 for upper bound of resistivity
## no bound case
ibound             !0
## lower bound case
##iboud              !1
##lower bo.log10[O.m]!1.0
## upper bound case
##iboud              !2
##upper bo.log10[O.m]!2.8
## nslice: number of slices you want to generate
nslice             !2
##---- Details for the first slice ------
## First slice file will be ${ID}_01.dat
## slice y = 0 (WE vertical slice)
A                  !$A
B                  !$B
D                  !$D
boundary left      !$bl
boundary right     !$br
boundary bottom    !$bb
boundary top       !$bt
##--- Details for the second slice -----
## second slice file will be ${ID}_02.dat
## slice x = 0 (SN vertical slice)
A                  !$A2
B                  !$B2
D                  !$D2
boundary left      !$bl
boundary right     !$br
boundary bottom    !$bb
boundary top       !$bt
##------------------ control file for slice.exe ends here !! ----------------------
EOF

#[2]## generate tmp.dat
if [ $igendat -eq 1 ];then
SRC="../src"
cd $SRC
#make clean
make
cd -
${SRC}/slice.exe << EOF
tmp.ctl
EOF
fi

#[3]## Draw picture
scl=8/5
range=$bl/$br/$bb/$bt
range2=0/8/0/5
scl2=20/20
range3=0/20/0/20
bb1=a0.25/a0.25Wesn
bb2=a1f0.25/a0.25WeSn
echo range = $range

gmt begin slice pdf

#gmt makecpt -Crainbow -T1.7/2.2/0.01 -I
gmt makecpt -Crainbow -T1.0/2.2/0.01 -I
#gmt makecpt -Crainbow -T0/3.0/0.01 -I

# WE slice
gmt basemap -JX$scl -R$range -Bxa0.25We+l"Easting [km]" -Bya0.25Sn+l"Elevation [km]" -Y8
gmt plot "${ID}_01.dat" -C -L
gmt plot -W0.8,black,- -L << EOF  # outline of G3 anomaly 
0.15  0.50
0.15  1.05
-0.10 1.05
-0.10 1.15
-0.20 1.15
-0.20 0.95
-0.15 0.95
-0.15 0.50
EOF
echo "0.5  4.5 $ID WE" | gmt text -JX$scl -R$range2 -G255 -F+f10p,Helvetica,black+jLM

# SN slice
gmt basemap -JX$scl -R$range -Bxa0.25We+l"Northing [km]" -Bya0.25Sn+l"Elevation[ km]" -Y-6.5
gmt plot "${ID}_02.dat" -C -L
gmt psxy  -W0.8,black,- -L << EOF # G3 anomaly
0.15  0.50
0.15  1.05
-0.15 1.05
-0.15 0.50
EOF
echo "0.5  4.5 $ID SN" | gmt text -JX$scl -R$range2  -Gwhite -F+f10p,Helvetica,black+jLM

gmt psscale -Dx10/0+w10/0.3 -B0.1  -C -X-1 #-G-0.3/2.5

gmt text -R0/20/0/20 -JX20/20 -F+f13p,Helvetica,black+jCM << EOF 
10.5 11.5 log@-10@- @%30%@~r@~@%%
10.5 10.8 (@~W@~m)
EOF

gmt end

open slice.pdf

