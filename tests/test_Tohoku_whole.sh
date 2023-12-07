#! /bin/sh
# test Tohoku

source $HOME/jenkins/TMTGEM/tests/func.sh

testMktopoTohoku(){
   mktopo Tohoku
   assertEquals 1 $?
}

testMeshgenTohoku(){
   meshgen Tohoku
   assertEquals 1 $?
}

testComcotTohoku(){
   runcomcot Tohoku
   assertEquals 1 $?
}

testTmtgemTohoku(){
   emrun Tohoku em b14
   assertEquals 1 $?
}

# load shunit2
. /usr/local/bin/shunit2