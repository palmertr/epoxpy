#!/bin/bash
#SBATCH -p gpuq 
#SBATCH --job-name=epoxpy
#SBATCH -n 14
#SBATCH --output=out.o
#SBATCH --mail-type=All
#SBATCH --mail-user=stephenthomas1@boisestate.edu
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1
module purge
#module use /scratch/erjank_project/mike_modules/modulefiles/
#module load hoomd/2.1.0-sp
module load cuda80/toolkit/8.0.44 
source activate dybond
export PYTHONPATH=$PYTHONPATH:/home/sthomas/scratch/projects/hoomd-blue/build

python run.py
