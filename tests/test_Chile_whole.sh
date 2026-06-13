#! /bin/sh
# test Chile

source tests/func.sh

testMktopoChile(){
   mktopo Chile
   assertEquals 0 $?
}

testMeshgenChile(){
   meshgen Chile
   assertEquals 0 $?
}

testComcotChile(){
   runcomcot Chile
   assertEquals 0 $?
}

testTmtgemChile(){
   emrun Chile em A01
   assertEquals 0 $?
}

# load shunit2
. /usr/local/bin/shunit2