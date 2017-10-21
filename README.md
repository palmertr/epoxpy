[![DOI](https://zenodo.org/badge/106355819.svg)](https://zenodo.org/badge/latestdoi/106355819)
[![Conda Badge](https://anaconda.org/cmelab/epoxpy/badges/version.svg)](https://anaconda.org/cmelab/epoxpy)

# About #
epoxpy is a python package for simulating polymerization in epoxy blends. It uses the [dybond plugin](https://bitbucket.org/cmelab/hoomd_blue) in HOOMD-blue to achieve this. We recommend installing using conda.

# Dependencies
* Required
	* Python >= 2.7
	* Numpy  >= 1.13.3
	* HOOMD-Blue ([dynamic_bonding](https://bitbucket.org/cmelab/hoomd_blue) branch)
	* [mbuild](http://mosdef-hub.github.io/mbuild/) == 0.6.1 
* Optional
	* [signac](http://signac.readthedocs.io/en/latest/)
	* [signac-flow](https://signac-flow.readthedocs.io/en/latest/)
	* [freud](http://glotzerlab.engin.umich.edu/freud/)
	
# Installation

epoxpy can be easily installed through conda. Follow instructions [here](https://conda.io/docs/install/quick.html#miniconda-quick-install-requirements) to install miniconda.

### Install using conda

Conda version 4.3 or greater is supported.
```
conda create --name epoxpy python=3.5
source activate epoxpy
conda install -c cmelab -c glotzer -c mosdef -c omnia epoxpy
```
#### Sample script
```
import epoxpy.abc_type_epoxy_simulation as es
import epoxpy.temperature_profile_builder as tpb

# define mixing parameters
mix_time = 1000
mix_kt = 2.0
# define curing parameters
cure_kt = 2.0
cure_time = 5e4
# define cure temperature profile
temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
temp_profile.add_state_point(cure_time, cure_kt)
# create an epoxy simulation object of type "ABC" which consists of 'A'
# particles, 'B' particles and 'C' particles mixed at a default ratio of 1:2:2
myEpoxySim = es.ABCTypeEpoxySimulation(sim_name='test',
                                       mix_time=mix_time,
                                       mix_kt=mix_kt,
                                       temp_prof=temp_profile,
                                       n_mul=2.0)# n_mul of 1 is 50 particles
myEpoxySim.execute()
```
#### Running jobs on clusters (or locally)

It is convenient to use the signac-flow based project called [epoxpy-flow](https://bitbucket.org/cmelab/epoxpy-flow) to run epoxpy jobs. 


### Install from source

```
git clone git@bitbucket.org:cmelab/epoxpy.git
cd epoxy_sim
conda env create -f conda_env.yml
pip install .
```

To check if the install from source was successful, run the tests. Do this also when you make changes and want to submit a pull request.

#### Running tests
```
cd epoxpy
pytest
```

#Contributing

To contribute, send a pull request to the "dev" branch as shown below:
```
git checkout dev
git checkout -b new_feature
# make the code changes and commit new branch
```
Now, create a pull request to the "dev" branch. This can be done in the bitbucket website.

# Maintainers

* Stephen Thomas (stephenthomas1@boisestate.edu)
* Mike Henry (mikehenry@boisestate.edu)
* Monet Alberts (monetalberts@u.boisestate.edu)
