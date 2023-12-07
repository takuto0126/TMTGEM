#! /bin/sh
# test Easter

source $HOME/jenkins/TMTGEM/tests/func.sh

testMktopoEaster(){
   mktopo Easter
   assertEquals 1 $?
}

testMeshgenEaster(){
   meshgen Easter
   assertEquals 1 $?
}

testComcotEaster(){
   runcomcot Easter
   assertEquals 1 $?
}

testTmtgemEaster(){
   emrun Easter em IPM
   assertEquals 1 $?
}

# load shunit2
. /usr/local/bin/shunit2