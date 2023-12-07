#! /bin/sh
# test Tohoku

TEST_FLDR=`pwd`

source $HOME/jenkins/TMTGEM/tests/func.sh
setUP(){
   echo "setUP is called"
}
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

tearDown(){
   echo "tearDown is called"
   cd $TEST_FLDR
}

# load shunit2
. /usr/local/bin/shunit2