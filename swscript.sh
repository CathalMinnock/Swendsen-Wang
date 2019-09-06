#!/bin/sh
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH -p compute
#SBATCH -J FINALq3s16n8

n=8
samples=10000
skip_steps=50
thermalize_steps=0
q=3
size=16
temps=40
beta_start=0
beta_end=1

beta_diff=`echo "scale=4; $beta_end-$beta_start" | bc`
beta_increment=`echo "scale=4; $beta_diff/$temps" | bc`
beta=$beta_start
directory_name="SW_q_${q}_size_${size}_n_${n}_therm_${thermalize_steps}"
rm -r $directory_name
mkdir $directory_name
for ((i = 0; i <= $temps; i++ ))
  do
    echo "HELLO" 
    filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    mpirun -n $n ./swprog -x $size -y $size -z $size -q $q -b $beta -f $filename -s $samples -a $skip_steps -t $thermalize_steps
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
for ((i = 0; i <= $temps; i++ ))
  do 
    input_filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    output_filename="${directory_name}/${directory_name}_beta_${beta}_autocorrelation.txt"
    rm $output_filename
    for((j = 0; j <= 30; j++ ))
      do
    	python autocorrelation.py $input_filename $output_filename $j
      done
    beta=`echo "scale=4; $beta+$beta_increment" | bc`
 done
 
 
