#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=r_output_%J_%a.txt
#SBATCH --error=r_error_%J_%a.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=cbcl_julia_submit
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/hurdleModelExplore
#SBATCH --array=1,2,3,4,5,6,7,8,9

module load Julia/1.7.2

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"


julia scripts/juliaCode/mHurdleFlex.jl data/data/CBCL_scale_${SLURM_ARRAY_TASK_ID}_Resp.csv data/data/CBCL_scale_${SLURM_ARRAY_TASK_ID}_Tabs.csv 