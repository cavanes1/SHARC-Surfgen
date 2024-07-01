#!/bin/bash
  
i=1
m=10

rm output.txt

while [ $i -le $m ]
do
 cd traj$i
 echo "entering traj$i"
 echo "entering traj$i" >> ../output.txt
 cp ../analyzegeom.py .
 python analyzegeom.py >> ../output.txt
 cd ..
 ((i++))
done

python readcollect.py

echo work complete
