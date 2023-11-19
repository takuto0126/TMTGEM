# coded on 2018.11.13
#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

ctl=comcot.ctl
src="../../comcotv1_7"

cd $src
make clean
make
cd -

time ${src}/comcot.exe < ${ctl}


