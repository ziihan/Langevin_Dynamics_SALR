#!/bin/bash
ND=256
Nthread=16
epsilon=6.0
for i in {1..1}
do
nohup ./BDkickoff_$ND-$i $i $Nthread $epsilon  > ./JobReport/BDkickoff_$ND-$i-$Nthread.opt&
echo "JOB $i has been submitted."
sleep 1s
done

