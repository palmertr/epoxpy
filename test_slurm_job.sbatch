#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name=epoxy_sim
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --output=log/main.o
#SBATCH --mail-type=All
#SBATCH --mail-user=stephenthomas1@boisestate.edu
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --gres=gpu:2

module purge
#module use /scratch/erjank_project/mike_modules/modulefiles/
#module load hoomd/2.1.0-sp
source activate epoxpy

python test_simulation.py
