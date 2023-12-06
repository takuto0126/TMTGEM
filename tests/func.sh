#!/bin/bash

function mktopo()
{
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/$1/topo/
  xyzfile=${fldr}/$2
  cd $fldr
  rm *.xyz #> /dev/null 2>&1
  chmod +x mktopo.sh
  ./mktopo.sh
  cd -
  if [ -e $xyzfile ];then
    return 1 
  else
    return 0 
  fi
}

function meshgen(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/$1/mesh/
  cd $fldr
  chmod +x clean.sh
  ./clean.sh
  chmod +x tetmeshgen.sh
  ./tetmeshgen.sh > log.txt
  cd -
  if [ -e ${fldr}em3d.msh ];then
    return 1 
  else
    return 0 
  fi
}

#-------------------------------------------- COMCOT
function sample3(){
  export PATH=$PATH:/usr/local/bin
  fldr=$HOME/jenkins/TMTGEM/Tohoku/flow/
  cd $fldr
  chmod +x clean.sh
  ./clean.sh
  chmod +x run_comcot.sh
  ./run_comcot.sh
  cd -
  if [ -e ${fldr}z_01_000300.dat ];then
    return 1 
  else
    return 0 
  fi
}

#------------------------------------------ em/run.sh
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
   echo $rms
   if [ `echo "$rms < 0.01" | bc` -eq 1 ]; then
    return 1 
  else
    return 0 
  fi
}