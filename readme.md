# epoxpy #


### What is this repository for? ###

* Used to run epoxy simulations using HOOMD as it's MD engine.
* Version 0.1

### How do I get set up? ###

* Install miniconda using instructions [here](https://conda.io/docs/install/quick.html#miniconda-quick-install-requirements).

* Clone the epoxpy repository and navigate into the folder
```
git clone git@bitbucket.org:cmelab/epoxy_sim.git
cd epoxy_sim
```

* Create the conda environment from the file 'conda_env.yml'
```
conda env create -f conda_env.yml
```
* The condo create should take care of all the dependencies except installing the epoxpy package. This is done by typing the following line in the terminal in the "epoxy_sim" directory.
```
pip install . --user
```
Alternatively, if you just want to install epoxpy package without cloning the repository, do
```
pip install -e git+ssh://git@bitbucket.org/cmelab/epoxy_sim.git@SimlationClass#egg=epoxpy --user
```

* How to run tests
```
python test_simulation.py
```
If you need to run the test in a different location like the scratch folder using slurm, copy test_simulation.py and test_slurm_job.sbatch into the desired folder and do
```
sbatch test_slurm_job.sbatch
```


### Who do I talk to? ###

* Repo owner or admin
* Other team members: stephenthomas1@boisestate.edu