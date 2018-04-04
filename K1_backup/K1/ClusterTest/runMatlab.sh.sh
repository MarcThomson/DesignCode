#!/bin/bash
#SBATCH -J test1.log
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -o test1.log-%j.out
#SBATCH -e test1.log-%j.err
#SBATCH --qos=normal 

module load matlab/2017b
matlab test1.m
