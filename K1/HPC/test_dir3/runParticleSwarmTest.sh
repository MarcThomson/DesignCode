#!/bin/bash
#SBATCH -J ParticleSwarmTest.log
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -o ParticleSwarmTest.log-%j.out
#SBATCH -e ParticleSwarmTest.log-%j.err
#SBATCH --qos=normal

cd /scratch/summit/math3916/Minimization/ParticleSwarm/
module load matlab
matlab -nosplash -nodesktop -r "clear; optimizationFile_ParticleSwarmTest;"
