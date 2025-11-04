#!/bin/bash
#
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=./outputText/r_output_%J_%a.txt
#SBATCH --error=./errorText/r_error_%J_%a.txt
#SBATCH --time=2:30:00
#SBATCH --job-name=hurdle_sim_submit
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/hurdleModelExplore
#SBATCH --array=1-1600

module purge
# module load GCC/11.3.0
# module load OpenMPI/4.1.4
# module load R/4.2.1
module load r/4.4.0

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
date
start=`date +%s`
Rscript ./scripts/rCode/informationUniqueExamine.R ${SLURM_ARRAY_TASK_ID}
echo "Done"
end=`date +%s`
runtime=$((end-start))
date
echo " Total runtime:" 
echo ${runtime}