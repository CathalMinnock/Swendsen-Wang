#!/bin/sh

n=4
q=2
size=16
temps=40
beta_start=0
beta_end=1
beta_increment=`echo "scale=4; $beta_end/$temps" | bc`
directory_name="SW_q_${q}_size_${size}_n_${n}"

beta=$beta_start
output_filename="${directory_name}_binder.txt"
rm $output_filename
for ((i = 0; i <= $temps; i++ ))
  do 
    input_filename="${directory_name}/${directory_name}_beta_${beta}.txt"
    python binder.py $input_filename $output_filename $beta
    beta=`echo "scale=4; $beta+$beta_increment" | bc`
 done
