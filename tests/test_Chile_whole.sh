#! /bin/sh
# test Chile

source tests/func.sh

testMktopoChile(){
   mktopo Chile
   assertEquals 1 $?
}

testMeshgenChile(){
   meshgen Chile
   assertEquals 1 $?
}

testComcotChile(){
   runcomcot Chile
   assertEquals 1 $?
}

testTmtgemChile(){
   emrun Chile em A01
   assertEquals 1 $?
}

# load shunit2
. /usr/local/bin/shunit2