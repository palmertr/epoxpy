# About #
epoxpy is a python package for simulating polymerization in epoxy blends. It uses the [dybond plugin](https://bitbucket.org/cmelab/hoomd_blue) in HOOMD-blue to achieve this.

### Prerequisites ###
* Required
	* Python >= 2.7
	* Numpy  >= 1.13.3
	* HOOMD-Blue ([dynamic_bonding](https://bitbucket.org/cmelab/hoomd_blue) branch)
	* [mbuild](http://mosdef-hub.github.io/mbuild/) == 0.6.1 
* Optional
	* [signac](http://signac.readthedocs.io/en/latest/)
	* [signac-flow](https://signac-flow.readthedocs.io/en/latest/)
	* [freud](http://glotzerlab.engin.umich.edu/freud/)
	
### Installation ###

epoxpy can be easily installed through [conda](https://conda.io/docs/install/quick.html#miniconda-quick-install-requirements). It is tested for python 3.5.

### Installing using conda

```
conda install -c cmelab epoxpy
```

### Installing from source

```
git clone git@bitbucket.org:cmelab/epoxpy.git
cd epoxy_sim
conda env create -f conda_env.yml
pip install .
```
### Running tests
```
cd epoxpy
pytest
```

### Running jobs on clusters (or locally)

It is convenient to use the signac-flow based project called [epoxpy-flow](https://bitbucket.org/cmelab/epoxpy-flow) to run epoxpy jobs. 

### Maintainers ###

* Stephen Thomas (stephenthomas1@boisestate.edu)
* Mike Henry (mikehenry@boisestate.edu)
* Monet Alberts (monetalberts@u.boisestate.edu)