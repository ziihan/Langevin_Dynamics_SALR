#!/bin/bash

##########mkdir OutSpherePut$Lx$Ly$Lz
jutil env activate -p jics30 -A jics30
####jutil env activate -p cjics22 -A jics22
module load intel-para/2019a-mt

for i in {11..20}
do
sbatch subjureca$i.sh
sleep 1s
done
