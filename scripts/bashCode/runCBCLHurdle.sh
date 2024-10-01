#!/bin/bash
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --chdir=/home/rosena/hurdleModelExplore/
#SBATCH --output=/home/rosena/CBCL_call.txt
#SBATCH --ntasks=12
#SBATCH --exclusive
#SBATCH --mem=36G

module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1
module load Julia

## Now write some output
echo "Processing started"
date

Rscript scripts/rCode/readCBCLDat.R

echo "Processing ended"
date