#!/bin/bash
#
#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=./outputText/r_output_%J_%a.txt
#SBATCH --error=./errorText/r_error_%J_%a.txt
#SBATCH --time=1:00:00
#SBATCH --job-name=hurdle_sim_submit
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/hurdleModelExplore
#SBATCH --array=1-64800

module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1
module load Julia/1.7.2

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
Rscript ./scripts/rCode/estSimHurdle.r ${SLURM_ARRAY_TASK_ID}