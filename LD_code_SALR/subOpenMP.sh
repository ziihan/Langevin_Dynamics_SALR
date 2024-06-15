#!/bin/bash -x
#SBATCH --nodes=1
#####SBATCH --ntasks=24
########SBATCH --ntasks-per-node=24
#SBATCH --ntasks-per-node=24
#######SBATCH --cpus-per-task=1
#SBATCH --output=./JobReport/testop24h.%j
#SBATCH --error=./JobReport/tester24h.%j
#SBATCH --time=24:00:00
#SBATCH --partition=batch

module load intel-para/2016b-mt

NC=100
Nthread=24
for i in {100..100}
do
#######srun --exclusive -n 1 ./LANE_$NC-$i $i $Nthread & 
srun ./DimColl_$NC-$i $i $Nthread &
echo "Job $i has been submitted."
####sleep 1s
done
wait
