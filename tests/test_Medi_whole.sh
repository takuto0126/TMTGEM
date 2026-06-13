#! /bin/sh
# test Mediterranean

source tests/func.sh

testMktopoMediterranean(){
   mktopo Mediterranean
   assertEquals 0 $?
}

testMeshgenMediterranean(){
   meshgen Mediterranean
   assertEquals 0 $?
}

testComcotMediterranean(){
   runcomcot Mediterranean
   assertEquals 0 $?
}

testTmtgemMediterranean(){
   emrun Mediterranean em A01
   assertEquals 0 $?
}

# load shunit2
. /usr/local/bin/shunit2