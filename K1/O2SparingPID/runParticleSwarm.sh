#!/bin/bash
#SBATCH -J ParticleSwarm.log
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -o ParticleSwarm.log-%j.out
#SBATCH -e ParticleSwarm.log-%j.err
#SBATCH --qos=normal

cd /scratch/summit/math3916/Minimization/ParticleSwarm/
module load matlab
matlab -nosplash -nodesktop -r "clear; optimizationFile_ParticleSwarm;"
