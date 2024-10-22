#!/bin/bash

## This script will take in one number, and run the next five increments of that number through the 
## mplus sim hurdle scirpt

int_num=${1}
end_val=$(( int_num + 2))
end_val_loop=$(( int_num + 1))
seq_vals=`seq ${int_num} ${end_val_loop}` 
echo ${seq_vals}
start=`date +%s`
for i in ${seq_vals}; do
  Rscript scripts/rCode/runMPlusHurdleSimSingle.R ${i} & 
done
Rscript scripts/rCode/runMPlusHurdleSimSingle.R ${end_val}
end=`date +%s`
runtime=$((end-start))
date
echo " Total runtime:" 
echo ${runtime}