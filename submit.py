import os
import sys
import shutil
import subprocess as sp

def slurm_job(email, job_time, sim_dir, queue="batch", job_name="epoxy_sim"): #, run_dir_0="foo_0/", run_dir_1="foo_1/"):
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
    job_string +="#SBATCH --gres=gpu:16\n"
    job_string +="\n"
    job_string +="module purge\n"
    # These next two lines need to be modified for your module env
    job_string +="module use /home/xsede/projects/p-dmr140097/mike/mod/modulefiles\n"
    job_string +="module load hoomd2.1.1-80-gc62a2fc\n"
    job_string +="cd {}\n".format(sim_dir)
    # This line assumes a 48hr wall clock time and may be commented out
    job_string +="export HOOMD_WALLTIME_STOP=$((`date +%s` + 48 * 3600 - 10 * 60))\n"
    job_string +="\n"
    job_string +="mpirun -np 1 --bind-to core --cpu-set 0  python {}sim.py {} {} --gpu=0  > {}{}job_0.o &\n".format(run_dir_0,  run_dir_0,  sys.argv[1], sim_dir, run_dir_0)
    job_string +="mpirun -np 1 --bind-to core --cpu-set 1  python {}sim.py {} {} --gpu=1  > {}{}job_1.o &\n".format(run_dir_1,  run_dir_1,  sys.argv[2], sim_dir, run_dir_1)
    job_string +="mpirun -np 1 --bind-to core --cpu-set 2  python {}sim.py {} {} --gpu=2  > {}{}job_0.o &\n".format(run_dir_2,  run_dir_2,  sys.argv[3], sim_dir, run_dir_2)   
    job_string +="mpirun -np 1 --bind-to core --cpu-set 3  python {}sim.py {} {} --gpu=3  > {}{}job_1.o &\n".format(run_dir_3,  run_dir_3,  sys.argv[4], sim_dir, run_dir_3)   
    job_string +="mpirun -np 1 --bind-to core --cpu-set 4  python {}sim.py {} {} --gpu=4  > {}{}job_0.o &\n".format(run_dir_4,  run_dir_4,  sys.argv[5], sim_dir, run_dir_4)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 5  python {}sim.py {} {} --gpu=5  > {}{}job_1.o &\n".format(run_dir_5,  run_dir_5,  sys.argv[6], sim_dir, run_dir_5)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 6  python {}sim.py {} {} --gpu=6  > {}{}job_0.o &\n".format(run_dir_6,  run_dir_6,  sys.argv[7], sim_dir, run_dir_6)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 7  python {}sim.py {} {} --gpu=7  > {}{}job_1.o &\n".format(run_dir_7,  run_dir_7,  sys.argv[8], sim_dir, run_dir_7)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 8  python {}sim.py {} {} --gpu=8  > {}{}job_0.o &\n".format(run_dir_8,  run_dir_8,  sys.argv[9], sim_dir, run_dir_8)          
    job_string +="mpirun -np 1 --bind-to core --cpu-set 9  python {}sim.py {} {} --gpu=9  > {}{}job_1.o &\n".format(run_dir_9,  run_dir_9,  sys.argv[10], sim_dir, run_dir_9)
    job_string +="mpirun -np 1 --bind-to core --cpu-set 10 python {}sim.py {} {} --gpu=10 > {}{}job_0.o &\n".format(run_dir_10, run_dir_10, sys.argv[11], sim_dir, run_dir_10)   
    job_string +="mpirun -np 1 --bind-to core --cpu-set 11 python {}sim.py {} {} --gpu=11 > {}{}job_1.o &\n".format(run_dir_11, run_dir_11, sys.argv[12], sim_dir, run_dir_11)   
    job_string +="mpirun -np 1 --bind-to core --cpu-set 12 python {}sim.py {} {} --gpu=12 > {}{}job_0.o &\n".format(run_dir_12, run_dir_12, sys.argv[13], sim_dir, run_dir_12)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 13 python {}sim.py {} {} --gpu=13 > {}{}job_1.o &\n".format(run_dir_13, run_dir_13, sys.argv[14], sim_dir, run_dir_13)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 14 python {}sim.py {} {} --gpu=14 > {}{}job_0.o &\n".format(run_dir_14, run_dir_14, sys.argv[15], sim_dir, run_dir_14)      
    job_string +="mpirun -np 1 --bind-to core --cpu-set 15 python {}sim.py {} {} --gpu=15 > {}{}job_1.o &\n".format(run_dir_15, run_dir_15, sys.argv[16], sim_dir, run_dir_15)      
    
    job_string +="wait\n"
    job_string +="echo 'all done!'\n"
    job_string +="exit"
    return job_string


def write_job_string(job_string, run_dir):
    with open(run_dir + "submit.sh", 'w') as out:
            out.write(job_string)


if __name__ == "__main__":
    email = "mikehenry@boisestate.edu"
    job_time = "00:30:00"
    # This should be the folder that sim.py, init.py, and submit.py are in
    sim_dir = "/cstor/xsede/projects/p-dmr140097/mike_epoxy_sim/"
    project_name = "msd-bonding-1e3"
    # This will be a sub folder in the sim_dir directory
    run_dir_0 = "runs/{}-kT-{}/".format(project_name, sys.argv[1])
    run_dir_1 = "runs/{}-kT-{}/".format(project_name, sys.argv[2])
    run_dir_2 = "runs/{}-kT-{}/".format(project_name, sys.argv[3])
    run_dir_3 = "runs/{}-kT-{}/".format(project_name, sys.argv[4])
    run_dir_4 = "runs/{}-kT-{}/".format(project_name, sys.argv[5])
    run_dir_5 = "runs/{}-kT-{}/".format(project_name, sys.argv[6])
    run_dir_6 = "runs/{}-kT-{}/".format(project_name, sys.argv[7])
    run_dir_7 = "runs/{}-kT-{}/".format(project_name, sys.argv[8])
    run_dir_8 = "runs/{}-kT-{}/".format(project_name, sys.argv[9])
    run_dir_9 = "runs/{}-kT-{}/".format(project_name, sys.argv[10])
    run_dir_10 = "runs/{}-kT-{}/".format(project_name, sys.argv[11])
    run_dir_11 = "runs/{}-kT-{}/".format(project_name, sys.argv[12])
    run_dir_12 = "runs/{}-kT-{}/".format(project_name, sys.argv[13])
    run_dir_13 = "runs/{}-kT-{}/".format(project_name, sys.argv[14])
    run_dir_14 = "runs/{}-kT-{}/".format(project_name, sys.argv[15])
    run_dir_15 = "runs/{}-kT-{}/".format(project_name, sys.argv[16])
    #make some dirs
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_0)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_1)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_2)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_3)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_4)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_5)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_6)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_7)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_8)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_9)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_10)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_11)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_12)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_13)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_14)
    sp.call(cmd, shell=True)
    cmd = "mkdir -p {}{}".format(sim_dir, run_dir_15)
    sp.call(cmd, shell=True)

    cwd = os.getcwd()
    cwd+= "/"
    job_string = slurm_job(email, job_time, sim_dir, queue="normal", job_name="epoxy_sim")#, run_dir_0=run_dir_0, run_dir_1 =  run_dir_1)
    
    write_job_string(job_string, sim_dir + run_dir_0)
    write_job_string(job_string, sim_dir + run_dir_1)
    write_job_string(job_string, sim_dir + run_dir_2)
    write_job_string(job_string, sim_dir + run_dir_3)
    write_job_string(job_string, sim_dir + run_dir_4)
    write_job_string(job_string, sim_dir + run_dir_5)
    write_job_string(job_string, sim_dir + run_dir_6)
    write_job_string(job_string, sim_dir + run_dir_7)
    write_job_string(job_string, sim_dir + run_dir_8)
    write_job_string(job_string, sim_dir + run_dir_9)
    write_job_string(job_string, sim_dir + run_dir_10)
    write_job_string(job_string, sim_dir + run_dir_11)
    write_job_string(job_string, sim_dir + run_dir_12)
    write_job_string(job_string, sim_dir + run_dir_13)
    write_job_string(job_string, sim_dir + run_dir_14)
    write_job_string(job_string, sim_dir + run_dir_15)

    shutil.copy("sim.py", sim_dir + run_dir_0 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_1 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_2 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_3 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_4 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_5 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_6 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_7 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_8 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_9 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_10 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_11 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_12 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_13 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_14 + "sim.py")
    shutil.copy("sim.py", sim_dir + run_dir_15 + "sim.py")

    shutil.copy("init.py", sim_dir + run_dir_0 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_1 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_2 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_3 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_4 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_5 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_6 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_7 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_8 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_9 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_10 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_11 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_12 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_13 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_14 + "init.py")
    shutil.copy("init.py", sim_dir + run_dir_15 + "init.py")

    cmd = "sbatch " + sim_dir + run_dir_0 + "submit.sh"
    sp.call(cmd, shell=True)
