import os
from simulation import Simulation
from bonding import Bonding
import hoomd
from hoomd import md
from hoomd import deprecated
from hoomd import dump
from abc import ABCMeta, abstractmethod


class EpoxySimulation(Simulation):
    """Common base class for all epoxy simulations.
          sim_name : name of simulation
          mix_time : number of time steps to run the equilibration (NVE simulation to get a nicely "shaken" phase)
          md_time  : number of time steps to run the molecular dynamics simulation (NVE)
          mix_kt   : temperature at which the mixing should be performed (during this time, we do an NVT)
          temp_prof: Dictionary of temperature and time which specifies the temperature profile to maintain during the
                     md_time. If md_time exceeds the last time step mentioned in the temp_prof, the last temperature is
                     maintained. If md_time is lesser than the last time step in the profile, the simulation ends without
                     completing the prescribed profile.
          log_write: time interval with which to write the log file
          dcd_write: time interval with which to write the dcd file
          output_dir: default is the working directory
          bond     : boolean value denoting whether to run the bonding routine for A's and B's
          bond_period: time interval between calls to the bonding routine
    """
    __metaclass__ = ABCMeta
    engine_name = 'HOOMD'

    def __init__(self, sim_name, mix_time, mix_kt, temp_prof, log_write=100, dcd_write=100, output_dir=os.getcwd(),
                 bond=False, bond_period=1e1, box=[3, 3, 3], dt=1e-2, legacy_bonding=False):
        Simulation.__init__(self, self.engine_name)
        self.simulation_name = sim_name
        self.mix_time = mix_time
        self.output_dir = os.path.join(output_dir,sim_name)
        self.bond = bond
        self.bond_period = bond_period
        self.mix_kT = mix_kt
        final_time = temp_prof.get_total_sim_time()
        md__total_time = final_time - mix_time
        self.md_time = md__total_time
        print('md time: {}'.format(self.md_time))
        self.temp_prof = temp_prof
        self.log_write = log_write
        self.dcd_write = dcd_write
        self.box = box
        self.msd_groups = None
        self.system = None
        self.dt = dt
        self.legacy_bonding = legacy_bonding

    def get_sim_name(self):
        return self.simulation_name

    def run_initial_structure(self):
        md.integrate.mode_standard(dt=self.dt)
        md.integrate.nve(group=hoomd.group.all())
        hoomd.run(self.mix_time)
        deprecated.dump.xml(group=hoomd.group.all(),
                            filename=os.path.join(self.output_dir, 'mix.hoomdxml'), all=True)

    def run_md(self):
        hoomd.run(self.md_time)
        deprecated.dump.xml(group=hoomd.group.all(),
                            filename=os.path.join(self.output_dir,
                                                  'final.hoomdxml'), all=True)

    @abstractmethod
    def set_initial_structure(self):
        pass

    @abstractmethod
    def set_force_field(self):
        pass

    @abstractmethod
    def setup_md_run(self):
        pass

    def configure_outputs(self):
        print('Configuring outputs. output_dir: {}'.format(self.output_dir))
        print('log_write: {} dcd_write: {}'.format(self.log_write, self.dcd_write))
        hoomd.meta.dump_metadata(filename=os.path.join(self.output_dir,
                                                       'metadata.json'), indent=2)
        deprecated.dump.xml(group=hoomd.group.all(),
                            filename=os.path.join(self.output_dir,
                                                  'start.hoomdxml'), all=True)
        hoomd.analyze.log(filename=os.path.join(self.output_dir, 'out.log'),
                          quantities=["pair_dpd_energy", "volume", "momentum", "potential_energy", "kinetic_energy",
                                      "temperature", "pressure", "bond_harmonic_energy"], period=self.log_write,
                          header_prefix='#', overwrite=True)
        dump.dcd(filename=os.path.join(self.output_dir, 'traj.dcd'), period=self.dcd_write, overwrite=True)
        dump.gsd(filename=os.path.join(self.output_dir, 'data.gsd'), period=self.dcd_write, group=hoomd.group.all(),
                 overwrite=True, static=['attribute', 'topology'])

    def initialize(self):
        hoomd.context.initialize()
        print('Creating simulation folder: {}'.format(self.output_dir))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        print('Initializing {}'.format(self.simulation_name))
        self.set_initial_structure()
        self.set_force_field()
        self.configure_outputs()
        self.run_initial_structure()

    def run(self):
        self.setup_md_run()
        print('Running {}'.format(self.simulation_name))

        deprecated.analyze.msd(groups=self.msd_groups, period=self.log_write,
                               filename=os.path.join(self.output_dir,
                                                     'msd.log'), header_prefix='#')
        if self.bond is True:
            log = hoomd.analyze.log(filename=None, quantities=["temperature"], period=self.bond_period)
            bonding_callback = Bonding(system=self.system, epoxy_sim=self, log=log, legacy_bonding=self.legacy_bonding)
            bond_callback = hoomd.analyze.callback(callback=bonding_callback, period=self.bond_period)

        self.run_md()

        if self.bond is True:
            bond_callback.disable()
            log.disable()

    def output(self):
        import numpy as np
        data = np.genfromtxt(fname=os.path.join(self.output_dir, 'out.log'), skip_header=True)

    def execute(self):
        print('Executing {}'.format(self.simulation_name))
        self.initialize()
        self.run()
        self.output()
        print("Finished executing {}".format(self.simulation_name))
