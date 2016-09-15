import os
import subprocess as sp

def slurm_job(email, job_time, sim_dir, queue="batch", job_name="epoxy_sim"):
    # sim_dir needs to be the dir sim.py is in
    job_string = "#!/bin/bash -l\n"
    job_string +="#SBATCH -p {}\n".format(queue)
    job_string +="#SBATCH -J {}\n".format(job_name)
    job_string +="#SBATCH -o job.o\n"
    job_string +="#SBATCH -N 1\n"
    job_string +="#SBATCH -n 16\n"
    job_string +="#SBATCH --mail-type=All\n"
    job_string +="#SBATCH --mail-user={}\n".format(email)
    job_string +="#SBATCH -t {}\n".format(job_time)
    job_string +="#SBATCH --exclusive\n"
    job_string +="#SBATCH --gres=gpu:2\n"
    job_string +="\n"
    job_string +="module purge\n"
    job_string +="module use /scratch/erjank_project/mike_modules/modulefiles/\n"
    job_string +="module load hoomd/0e1d572 slurm\n"
    job_string +="cd {}\n".format(sim_dir)
    # This next line is not used yet
    job_string +="export HOOMD_WALLTIME_STOP=$((`date +%s` + 12 * 3600 - 10 * 60))\n"
    job_string +="\n"
    job_string +="mpirun -np 1 --bind-to core --cpu-set 0 python sim.py --gpu=0 > job_a.o &\n"
    job_string +="mpirun -np 1 --bind-to core --cpu-set 1 python sim.py --gpu=1 > job_b.o &\n"
    job_string +="wait"
    return job_string


def write_job_string(job_string):
    with open("submit.sh", 'w') as out:
            out.write(job_string)


if __name__ == "__main__":
    email = "mikehenry@boisestate.edu"
    job_time = "1:00:00"
    sim_dir = "/scratch/erjank_project/mike_epoxy_sim/"

    job_string = slurm_job(email, job_time, sim_dir, queue="batch", job_name="epoxy_sim")
    write_job_string(job_string)
    cmd = "sbatch submit.sh"
    sp.run(cmd, shell=True)
