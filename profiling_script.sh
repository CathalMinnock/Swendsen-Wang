#!/bin/sh
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH -p compute
#SBATCH -J FINALq3s16n8

n=8
samples=1000
size=8

skip_steps=1
q=3
beta=0.4

module load cports
module load scorep
export SCOREP_ENABLE_TRACING=true
export SCOREP_ENABLE_PROFILING=true
export SCOREP_METRIC_PAPI=PAPI_L2_DCM
export SCOREP_TOTAL_MEMORY=3G
make
mpirun -n $n ./swprog -x $size -y $size -z $size -q $q -b $beta -f $filename -s $samples -a $skip_steps
