#! /bin/sh
# test Chile

function sample(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/Chile/topo/
  cd $fldr
  rm *.xyz #> /dev/null 2>&1
  chmod +x mktopo.sh
  ./mktopo.sh
  cd -
  if [ -e ${fldr}/W-81E-69S-53N-41_1min.xyz ];then
#  if [ 1 == 1 ];then
    echo 1 > result.txt
  else
    echo 0 > result.txt
  fi
}

testMktopoChile()
{
   sample
   o1=`head -1 result.txt`
   assertEquals 1 $o1
#   assertEquals 1 1
}

function sample2(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/Chile/mesh/
  cd $fldr
  chmod +x clean.sh
  ./clean.sh
  chmod +x tetmeshgen.sh
  ./tetmeshgen.sh > log.txt
  cd -
  if [ -e ${fldr}em3d.msh ];then
#  if [ 1 == 1 ];then
    echo 1 > result.txt
  else
    echo 0 > result.txt
  fi
}

testMeshgenChile()
{
   sample2
   o1=`head -1 result.txt`
   assertEquals 1 $o1
#   assertEquals 1 1
}

#-------------------------------------------- COMCOT
function sample3(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/Chile/flow/
  cd $fldr
  chmod +x clean.sh
  ./clean.sh
  chmod +x run_comcot.sh
  ./run_comcot.sh
  cd -
  if [ -e ${fldr}z_01_000300.dat ];then
#  if [ 1 == 1 ];then
    echo 1 > result.txt
  else
    echo 0 > result.txt
  fi
}

testComcotChile()
{
   sample3
   o1=`head -1 result.txt`
   assertEquals 1 $o1
#   assertEquals 1 1
}

function sample4(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/Chile/em/
  cd $fldr
  chmod +x clean.sh
  ./clean.sh
  chmod +x run.sh
  ./run.sh
  cd -
  fil1=~/jenkins/TMTGEM/Chile/em/bxyz/A01_bxyz_ts.dat
  fil2=Chile_ref/A01_bxyz_ts.dat
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

testTmtgemChile()
{
   sample4
   o1=`head -1 result.txt`
   assertEquals 1 $o1
#   assertEquals 1 1
}

# load shunit2
. /usr/local/bin/shunit2