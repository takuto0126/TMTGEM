#! /bin/sh
# test Mediterranean

source $HOME/jenkins/TMTGEM/tests/func.sh

testMktopoMediterranean(){
   mktopo Mediterranean
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
   emrun Mediterranean em A01
   assertEquals 1 $?
}

# load shunit2
. /usr/local/bin/shunit2