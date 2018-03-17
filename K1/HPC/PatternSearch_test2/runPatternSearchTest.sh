#!/bin/bash
#SBATCH -J PatternSearchTest.log
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -o PatternSearchTest.log-%j.out
#SBATCH -e PatternSearchTest.log-%j.err
#SBATCH --qos=normal

cd /scratch/summit/math3916/Minimization/PatternSearch/
module load matlab
matlab -nosplash -nodesktop -r "clear; optimizationFile_PatternSearchTest;"
