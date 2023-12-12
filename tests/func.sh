#!/bin/bash

TMTGEM_HOME=`pwd`
TEST_FLDR=${TMTGEM_HOME}/tests

setUp(){
  export PATH=$PATH:/usr/local/bin
  echo "setUp is called"
}

tearDown(){
   echo "tearDown is called"
   cd $TEST_FLDR
}


function mktopo()
{
  fldr=${TMTGEM_HOME}/$1/topo/
  xyzfile=${fldr}topo.xyz
  cd $fldr
  rm *.xyz #> /dev/null 2>&1
  chmod +x mktopo.sh
  ./mktopo.sh
#  cd -
  if [ -e $xyzfile ];then
    echo "mktopo SUCCESS!!"
    return 1 
  else
    echo "mktopo Failure..."
    return 0 
  fi
}

function meshgen(){
  export PATH=$PATH:/usr/local/bin
  fldr=${TMTGEM_HOME}/$1/mesh/
  cd $fldr
  chmod +x clean.sh
  chmod +x tetmeshgen.sh # clean.sh is included intetmeshgen.sh
  ./tetmeshgen.sh > log.txt
#  cd -
  if [ -e ${fldr}em3d.msh ];then
    echo "meshgen SUCCESS!!"
    return 1
  else
    echo "meshgen Failure..."
    return 0 
  fi
}

#-------------------------------------------- COMCOT
function runcomcot(){
  export PATH=$PATH:/usr/local/bin
  fldr=${TMTGEM_HOME}/$1/flow/
  cd $fldr
  chmod +x clean.sh
  chmod +x run_comcot.sh # clean.sh is included 2023.12.07
  ./run_comcot.sh
#  cd -
  if [ -e ${fldr}z_01_000300.dat ];then
    echo "runcomcot SUCCESS!!"
    return 1 
  else
    echo "runcomcot Failure..."
    return 0 
  fi
}

#------------------------------------------ em/run.sh
function emrun(){
  export PATH=$PATH:/usr/local/bin
  fldr=${TMTGEM_HOME}/$1/$2/
  compfile=$3_bxyz_ts.dat
  cd $fldr
  chmod +x clean.sh
  chmod +x run.sh # clean.sh is included in run.sh 2023.12.07
  ./run.sh
#  cd -
  fil1=${fldr}bxyz/$compfile
  fil2=${TMTGEM_HOME}/${1}_ref/$compfile
  rms=`paste $fil1 $fil2 | awk '{m+=($2 - $6)^2}END{printf "%15.7f", sqrt(m/NR);}'`
   echo RMS = $rms
   if [ `echo "$rms < 0.01" | bc` -eq 1 ]; then
    echo "emrun SUCCESS!!"
    return 1
  else
    echo "emrun Failure..."
    return 0 
  fi
}
