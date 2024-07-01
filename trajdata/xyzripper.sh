#/bin/bash

i=0
while [ $i -le 1000 ]
do
  let j=6*$i+1
  let k=6*$i+6
  sed -n "$j","$k"p output.xyz > $i.xyz
   ((i++))
done
