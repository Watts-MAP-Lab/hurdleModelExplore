#!/bin/bash

## This will be the wrapper for the runMplusSim.sh script
## It is going to try to make sure 5 instances of the mplus sim are running until completion

cd /home/arosen/Documents/hurdleModelExplore
all_int=`seq 1 3 100`
for l in ${all_int}; do
  /bin/bash scripts/bashCode/runMPlusSim.sh ${l}  > ~/Documents/hurdleModelExplore/${l}.txt ; 
done