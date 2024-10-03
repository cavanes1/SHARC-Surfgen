#/bin/bash

PRIMARYDIR=$(pwd)

cd QM
cd surface
#first read surfgen surface
#s0_s1 DPEM
cd s0_s1
cp ../../QM.in .
$PRIMARYDIR/../../sharc-surfgen/surfgen.x > log
cp e.out ../s0_s1.out
cd ..
#t1
cd t1
cp ../../QM.in .
$PRIMARYDIR/../../sharc-surfgen/surfgen.x > log
cp e.out ../t1.out
cd ..
#go into property directory
cd dip
cp ../../QM.in .
$PRIMARYDIR/../../sharc-surfgen/dip.x > log
cp dip.out ../dip.out
cd ..
#soc
cd soc
cp ../../QM.in .
$PRIMARYDIR/../../sharc-surfgen/qm.x > log
#cp QM.out to main QM dir
cp QM.out ../../.
cd ../../
