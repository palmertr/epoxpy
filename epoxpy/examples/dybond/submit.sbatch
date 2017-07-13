#!/bin/bash
#SBATCH -p gpuq 
#SBATCH --job-name=epoxpy
#SBATCH -n 14
#SBATCH --output=out.%jo
#SBATCH --mail-type=All
#SBATCH --mail-user=mikehenry@boisestate.edu
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1
module purge
module use /scratch/mhenry/mike_modules/modulefiles/
module load hoomd/dybond-hoomd
#module use /scratch/erjank_project/mike_modules/modulefiles/
#module load hoomd/2.1.0-sp
#module load cuda80/toolkit/8.0.44 
which python
export PATH="/home/$USER/miniconda3/bin:$PATH" # CHECK PATH PLZ
which python
source activate epoxpy_env
python -u run.py
