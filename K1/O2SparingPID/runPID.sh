#!/bin/bash
#SBATCH -J ParticleSwarm_PID.log
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -o ParticleSwarm_PID.log-%j.out
#SBATCH -e ParticleSwarm_PID.log-%j.err
#SBATCH --qos=normal

cd /scratch/summit/math3916/Minimization/PID
module load matlab
matlab -nosplash -nodesktop -r "clear; optimizationFile_PID;"
