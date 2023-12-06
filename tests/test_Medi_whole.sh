#! /bin/sh
# test Mediterranean

source $HOME/jenkins/TMTGEM/tests/func.sh

testMktopoMediterranean(){
   mktopo Mediterranean W130E155S33N45_1min.xyz
   assertEquals 1 $?
}

testMeshgenMediterranean(){
   meshgen Mediterranean
   assertEquals 1 $?
}

testComcotMediterranean(){
   runcomcot Mediterranean
   assertEquals 1 $?
}

testTmtgemMediterranean(){
   emrun Mediterranean em b14
   assertEquals 1 $?
}

# load shunit2
. /usr/local/bin/shunit2