#!/bin/sh
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -p compute
#SBATCH -J q3s16n16

n=16
samples=10000
q=3
size=16
temps=40
beta_start=0
beta_end=1
beta_increment=`echo "scale=4; $beta_end/$temps" | bc`
beta=$beta_start

directory_name="SW_q_${q}_size_${size}_n_${n}"
rm -r $directory_name
mkdir $directory_name
module load cports
module load openmpi
for ((i = 0; i <= $temps; i++ ))
  do
    echo "HELLO" 
    filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    mpirun -n $n ./swprog -x $size -y $size -z $size -q $q -b $beta -f $filename -s $samples
    beta=`echo "scale=4; $beta+$beta_increment" | bc`
 done
 
beta=$beta_start
output_filename="${directory_name}_magnetization.txt"
rm $output_filename
for ((i = 0; i <= $temps; i++ ))
  do 
    input_filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    python jackknife.py $input_filename $output_filename $beta
    beta=`echo "scale=4; $beta+$beta_increment" | bc`
 done
 
beta=$beta_start
gcc -Wall -o autocorrelation_prog autocorrelation.c
for ((i = 0; i <= $temps; i++ ))
  do 
    input_filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    output_filename="${directory_name}/${directory_name}_beta_${beta}_autocorrelation.txt"
    rm $output_filename
    for((j = 0; j <= 30; j++ ))
      do
    	./autocorrelation_prog -i $input_filename -o $output_filename -t $j
      done
    beta=`echo "scale=4; $beta+$beta_increment" | bc`
 done
 
 
