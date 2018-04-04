#!/bin/bash
#SBATCH -J PatternSearch.log
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -o PatternSearch.log-%j.out
#SBATCH -e PatternSearch.log-%j.err
#SBATCH --qos=normal

cd /scratch/summit/math3916/Minimization/PatternSearch/
module load matlab
matlab -nosplash -nodesktop -r "clear; optimizationFile_PatternSearch;"
