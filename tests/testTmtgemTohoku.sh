#! /bin/sh
# test Tohoku

#-------------------------------------------- Tohoku/EM
function sample4(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/Tohoku/em/
  cd $fldr
  chmod +x clean.sh
  ./clean.sh
  chmod +x run.sh
  ./run.sh
  cd -
  fil1=~/jenkins/TMTGEM/Tohoku/em/bxyz/b14_bxyz_ts.dat
  fil2=Tohoku_ref/b14_bxyz_ts.dat
  rms=`paste $fil1 $fil2 | awk '{m+=($2 - $6)^2}END{printf "%15.7f", sqrt(m/NR);}'`
#  if [ -e ${fldr}bxyz/bxyz_xy2D000060.dat ];then
#  if [ 1 == 1 ];then
   echo $rms
   if [ $rms < 0.01 ]; then
    echo 1 > result.txt
  else
    echo 0 > result.txt
  fi
}

testTmtgemTohoku()
{
   sample4
   o1=`head -1 result.txt`
   assertEquals 1 $o1
#   assertEquals 1 1
}

# load shunit2
. /usr/local/bin/shunit2