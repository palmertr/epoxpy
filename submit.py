import os
import sys
import shutil
import subprocess as sp

def slurm_job(email, job_time, sim_dir, queue="batch", job_name="epoxy_sim", run_dir_0="foo_0/", run_dir_1="foo_1/"):
    # sim_dir needs to be the dir sim.py is in
    job_string = "#!/bin/bash -l\n"
    if queue == "quick":
        job_string +="#SBATCH -A quick\n"
        job_time = "1:00:00"
    job_string +="#SBATCH -p {}\n".format(queue)
    job_string +="#SBATCH -J {}\n".format(job_name)
    job_string +="#SBATCH -N 1\n"
    job_string +="#SBATCH -n 16\n"
    job_string +="#SBATCH -o {}{}main.o\n".format(sim_dir, run_dir_0)
    job_string +="#SBATCH --mail-type=All\n"
    job_string +="#SBATCH --mail-user={}\n".format(email)
    job_string +="#SBATCH -t {}\n".format(job_time)
    job_string +="#SBATCH --exclusive\n"
    job_string +="#SBATCH --gres=gpu:2\n"
    job_string +="\n"
    job_string +="module purge\n"
    # These next two lines need to be modified for your module env
    job_string +="module use /scratch/erjank_project/mike_modules/modulefiles/\n"
    job_string +="module load hoomd/2.1.0\n"
    job_string +="cd {}\n".format(sim_dir)
    # This line assumes a 48hr wall clock time and may be commented out
    job_string +="export HOOMD_WALLTIME_STOP=$((`date +%s` + 48 * 3600 - 10 * 60))\n"
    job_string +="\n"
    job_string +="mpirun -np 1 --bind-to core --cpu-set 0 python {}sim.py {} {} --gpu=0 > {}{}job_0.o &\n".format(run_dir_0, run_dir_0, sys.argv[1], sim_dir, run_dir_0)
    job_string +="mpirun -np 1 --bind-to core --cpu-set 1 python {}sim.py {} {} --gpu=1 > {}{}job_1.o &\n".format(run_dir_1, run_dir_1, sys.argv[2], sim_dir, run_dir_1)
    job_string +="wait\n"
    job_string +="echo 'all done!'"
    return job_string


def write_job_string(job_string, run_dir):
    with open(run_dir + "submit.sh", 'w') as out:
            out.write(job_string)


if __name__ == "__main__":
    email = "mikehenry@boisestate.edu"
    job_time = "04:00:00"
    # This should be the folder that sim.py, init.py, and submit.py are in
    sim_dir = "/scratch/erjank_project/mike_epoxy_sim/"
    project_name = "test-msd-no-bonding"
    # This will be a sub folder in the sim_dir directory
    run_dir_0 = "runs/{}-kT-{}/".format(project_name, sys.argv[1])
    run_dir_1 = "runs/{}-kT-{}/".format(project_name, sys.argv[2])

    #make some dirs
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_0)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_1)
    sp.call(cmd, shell=True)

    cwd = os.getcwd()
    cwd+= "/"
    #run_name = "dpdc_debug_bonding_p10g0_{}/".format(run_name_postfix)
    job_string = slurm_job(email, job_time, sim_dir, queue="quick", job_name="epoxy_sim", run_dir_0=run_dir_0, run_dir_1 =  run_dir_1)
    write_job_string(job_string, run_dir_0)
    write_job_string(job_string, run_dir_1)

    shutil.copy("sim.py", cwd + run_dir_0 + "sim.py")
    shutil.copy("sim.py", cwd + run_dir_1 + "sim.py")

    shutil.copy("init.py", cwd + run_dir_0 + "init.py")
    shutil.copy("init.py", cwd + run_dir_1 + "init.py")


    cmd = "sbatch "+ run_dir_0 + "submit.sh"
    sp.call(cmd, shell=True)
