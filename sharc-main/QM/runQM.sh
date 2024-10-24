#/bin/bash

PRIMARYDIR=$(pwd)

cd QM
cd surface
#first read surfgen surface
cp ../QM.in .
$PRIMARYDIR/../../sharc-surfgen/surfgen.x > log
#cp QM.out to main QM dir
cp QM.out ../.
cd ..
