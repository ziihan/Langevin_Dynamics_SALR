#!/bin/bash
ND=256
mkdir Particlediffusion$ND
mkdir /p/scratch/jics30/BDfiles
mkdir /p/scratch/jics30/BDfiles/Particlediffusion$ND
mkdir /p/scratch/jics30/BDfiles/Particlediffusion$ND/salr_restart_files
module load intel-para/2019a-mt

for i in {16..20}
do
make clean
make -j  
cp bd_md LangevinQ2d_$ND-$i
####sleep 1s
done
