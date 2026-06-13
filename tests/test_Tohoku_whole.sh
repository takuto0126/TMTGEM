#! /bin/bash
# test Tohoku

source tests/func.sh

testMktopoTohoku(){
   mktopo Tohoku
   assertEquals 0 $?
}

testMeshgenTohoku(){
   meshgen Tohoku
   assertEquals 0 $?
}

testComcotTohoku(){
   runcomcot Tohoku
   assertEquals 0 $?
}

testTmtgemTohoku_em(){
   emrun Tohoku em b14
   assertEquals 0 $?
}

testTmtgemTohoku_IGRF(){
   emrun Tohoku em_IGRF b14
   assertEquals 0 $?
}

testTmtgemTohoku_woa(){
   emrun Tohoku em_woa b14
   assertEquals 0 $?
}


# load shunit2
. /usr/bin/shunit2