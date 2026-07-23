#!/bin/bash

TMTGEM_HOME=`pwd` # TMTGEM folder
TEST_FLDR=${TMTGEM_HOME}/tests

setUp(){
  export PATH=$PATH:/usr/local/bin
  #echo "setUp is called"
}

tearDown(){
   #echo "tearDown is called"
   cd $TEST_FLDR
}


function mktopo()
{
  fldr=${TMTGEM_HOME}/$1/topo/
  xyzfile=${fldr}topo.xyz
  logfile=${TMTGEM_HOME}/tests/log_${1}_mktopo.txt # 2026.06.14
  cd $fldr
  rm *.xyz 2>/dev/null
  chmod +x mktopo.sh

  echo "Running mktopo: $1"
  echo "Log: $logfile"
  ./mktopo.sh > $logfile 2>&1

#  cd -
  if [ -e $xyzfile ];then
    echo "mktopo SUCCESS!!"
    return 0 
  else
    echo "mktopo Failure..."
    echo "See log: $logfile"
    return 1
  fi
}

function meshgen(){
  export PATH=$PATH:/usr/local/bin
  fldr=${TMTGEM_HOME}/$1/mesh/
  logfile=${TMTGEM_HOME}/tests/log_${1}_meshgen.txt # 2026.06.14
  cd $fldr
  chmod +x clean.sh
  chmod +x tetmeshgen.sh # clean.sh is included intetmeshgen.sh
  ./tetmeshgen.sh > $logfile 2>&1
#  cd -
  if [ -e ${fldr}em3d.msh ];then
    echo "meshgen SUCCESS!!"
    return 0
  else
    echo "meshgen Failure..."
    echo "See log: $logfile"
    return 1 
  fi
}

#-------------------------------------------- COMCOT
function runcomcot(){
  export PATH=$PATH:/usr/local/bin
  fldr=${TMTGEM_HOME}/$1/flow/
  logfile=${TMTGEM_HOME}/tests/log_${1}_runcomcot.txt # 2026.06.14
  cd $fldr
  chmod +x clean.sh
  chmod +x run_comcot.sh # clean.sh is included 2023.12.07
  ./run_comcot.sh > $logfile 2>&1
#  cd -
  if [ -e ${fldr}z_01_000300.dat ];then
    echo "runcomcot SUCCESS!!"
    return 0 
  else
    echo "runcomcot Failure..."
    echo "See log: $logfile"
    return 1 
  fi
}

#------------------------------------------ em/run.sh
function emrun(){
  export PATH=$PATH:/usr/local/bin
  fldr=${TMTGEM_HOME}/$1/$2/
  compfile=$3_bxyz_ts.dat
  logfile=${TMTGEM_HOME}/tests/log_${1}_${2}_${3}.txt # 2026.06.14
  cd $fldr
  chmod +x clean.sh
  chmod +x run.sh # clean.sh is included in run.sh 2023.12.07
  echo "Running emrun: $1 $2 $3" # 2026.06.14
  echo "Log: $logfile"           # 2026.06.14
  ./run.sh > $logfile 2>&1 # 2026.06.14
  #  cd -
  fil1=${fldr}bxyz/$compfile
  fil2=${TMTGEM_HOME}/tests/${1}_ref/$compfile
  rms=`paste $fil1 $fil2 | awk '{m+=($2 - $6)^2}END{printf "%15.7f", sqrt(m/NR);}'`
   echo RMS = $rms
   if [ `echo "$rms < 0.01" | bc` -eq 1 ]; then
    echo "emrun SUCCESS!!"
    return 0
  else
    echo "emrun Failure..."
    echo "See log: $logfile"
    return 1 
  fi
}
