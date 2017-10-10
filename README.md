# epoxpy #
epoxpy is a python package for simulating polymerization in epoxy blends. It uses the [dybond plugin](https://bitbucket.org/cmelab/hoomd_blue) in HOOMD-blue to achieve this.

### Prerequisites ###
* Required
	* Python >= 2.7
	* Numpy  >= 1.13.3
	* HOOMD-Blue ([dynamic_bonding](https://bitbucket.org/cmelab/hoomd_blue) branch)
* Optional
	* [signac](http://signac.readthedocs.io/en/latest/)
	* [signac-flow](https://signac-flow.readthedocs.io/en/latest/)
	* [freud](http://glotzerlab.engin.umich.edu/freud/)
	* [mbuild](http://mosdef-hub.github.io/mbuild/)
	
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
cd epoxpy
pytest
```

* Running multiple jobs on clusters
It is convenient to use the signac-flow based project called [epoxpy-flow](https://bitbucket.org/cmelab/epoxpy-flow) to run epoxpy jobs. 

### Who do I talk to? ###

* Repo owner or admin
* Other team members: stephenthomas1@boisestate.edu