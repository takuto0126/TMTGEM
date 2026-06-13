#! /bin/sh
# test Easter

source tests/func.sh

testMktopoEaster(){
   mktopo Easter
   assertEquals 0 $?
}

testMeshgenEaster(){
   meshgen Easter
   assertEquals 0 $?
}

testComcotEaster(){
   runcomcot Easter
   assertEquals 0 $?
}

testTmtgemEaster(){
   emrun Easter em IPM
   assertEquals 0 $?
}

# load shunit2
. /usr/local/bin/shunit2