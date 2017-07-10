# epoxpy #


### What is this repository for? ###

* Used to run epoxy simulations using HOOMD as it's MD engine.
* Version 0.1

### How do I get set up? ###

*Install miniconda using instructions [here](https://conda.io/docs/install/quick.html#miniconda-quick-install-requirements).

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
pip install .
```
* How to run tests
```
python test_simulation.py

cd epoxpy
pytest
```

### Who do I talk to? ###

* Repo owner or admin
* Other team members: stephenthomas1@boisestate.edu

### Tips for Clusters ###

## R2 ##


```
#!bash

module use /scratch/mhenry/mike_modules/modulefiles/
module load hoomd/dybond-hoomd
```


## Kestrel ##


```
#!bash

module use /scratch/erjank_project/mike_modules/modulefiles/
module load hoomd/2.1.7-dybond
```