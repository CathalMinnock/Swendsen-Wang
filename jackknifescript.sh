#!/bin/sh
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -p compute
#SBATCH -J sw_q2

q=2
samples=1000
size=8
temps=20
beta_start=0
beta_end=1
beta_increment=`echo "scale=4; $beta_end/$temps" | bc`
beta=$beta_start

directory_name="SW_q_${q}_size_${size}"
output_filename="${directory_name}_magnetization.txt"
for ((i = 0; i <= $temps; i++ ))
  do 
  	input_filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    python jackknife.py $input_filename $output_filename $beta
    beta=`echo "scale=4; $beta+$beta_increment" | bc`
 done
