#!/bin/bash -l
#SBATCH -p quick
#SBATCH -J epoxy_sim
#SBATCH -o job.o
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=All
#SBATCH --mail-user=mikehenry@boisestate.edu
#SBATCH -t 1:00:00
#SBATCH --exclusive
#SBATCH --gres=gpu:2

module purge
module use /scratch/erjank_project/mike_modules/modulefiles/
module load hoomd/0e1d572 slurm
cd /scratch/erjank_project/mike_epoxy_sim/
export HOOMD_WALLTIME_STOP=$((`date +%s` + 12 * 3600 - 10 * 60))

mpirun -np 1 --bind-to core --cpu-set 0 python sim.py A --gpu=0 > job_a.o &
mpirun -np 1 --bind-to core --cpu-set 1 python sim.py B --gpu=1 > job_b.o &
wait