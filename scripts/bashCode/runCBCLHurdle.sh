#!/bin/bash
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --chdir=/home/rosena/hurdleModelExplore/
#SBATCH --output=/home/rosena/CBCL_call.txt
#SBATCH --ntasks=9
#SBATCH --exclusive
#SBATCH --mem=20G

module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1
module load Julia

## Now write some output
echo "Processing started"
date

#Rscript scripts/rCode/readCBCLDat.R

julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_1_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_1_Tabs.csv > ~/1.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_2_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_2_Tabs.csv > ~/2.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_3_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_3_Tabs.csv > ~/3.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_4_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_4_Tabs.csv > ~/4.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_5_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_5_Tabs.csv > ~/5.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_6_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_6_Tabs.csv > ~/6.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_7_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_7_Tabs.csv > ~/7.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_8_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_8_Tabs.csv > ~/8.txt &
julia /home/rosena/hurdleModelExplore/scripts/juliaCode /home/rosena/hurdleModelExplore/data/CBCL_scale_9_Resp.csv /home/rosena/hurdleModelExplore/data/CBCL_scale_9_Tabs.csv > ~/9.txt



echo "Processing ended"
date