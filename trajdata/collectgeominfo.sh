#!/bin/bash -l

#SBATCH --partition=shared
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --job-name=collectgeominfo
  
i=1
m=5000

rm output.txt
mkdir GEOMS

while [ $i -le $m ]
do
 cd traj$i
 echo "entering traj$i"
 echo "entering traj$i" >> ../output.txt
 cp ../analyzegeom.py .
 python analyzegeom.py $i 1 >> ../output.txt
 cp $i.xyz ../GEOMS
 cd ..
 ((i++))
done

python readcollect.py

echo work complete
