#/bin/bash

PRIMARYDIR=$(pwd)

cd QM
cd interface
cp ../QM.in .
$PRIMARYDIR/../../sharc-surfgen/sharc-surfgen.x > QM.log
cp QM.out ../.
cd ..
