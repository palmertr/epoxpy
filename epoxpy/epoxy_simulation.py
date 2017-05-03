import os
from epoxpy.simulation import Simulation
from epoxpy.bonding import LegacyBonding, FreudBonding
import hoomd
from hoomd import md
from hoomd import deprecated
from hoomd import dump
from abc import ABCMeta, abstractmethod
import random


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
          
          kwargs: These are parameters which might not be necessary from a user perspective. Includes parameters for 
          backward compatibility also.
                legacy_bonding: boolean value indicating whether to use legacy bonding or frued bonding routine
                                Default is "False"
                exclude_mixing_in_output: boolean value indicating whether the initial mixing phase should appear 
                                          in the output files. This may or may not be needed for analysis. 
                                          Default is "False"
                init_file_name: Full path to the initial file being written to disk during initialization. Please do
                                not use this to initialize from an initial structure created externally. That is not
                                the current intention of this argument. For now, use this to switch the use of gsd file
                                vs. hoomdxml file. Why did I not just use one of them instead? Because hoomdxml is
                                human readable and is handy for debugging, but is deprecated. So if they remove support
                                its easy to switch to using gsd file format.
                shrink: boolean value indicating whether the initial structure will be shrunk to reach the specified
                        density. default is "True"
                shrink_time: number of timesteps to run hoomd to shrink the initial volume to desired density.
                             Default is 1 time step.
    """
    __metaclass__ = ABCMeta
    engine_name = 'HOOMD'

    def __init__(self, sim_name, mix_time, mix_kt, temp_prof, log_write=100, dcd_write=100, output_dir=os.getcwd(),
                 bond=False, bond_period=1e1, box=[3, 3, 3], dt=1e-2, density=1.0, activation_energy=0.1):
        Simulation.__init__(self, self.engine_name)
        self.simulation_name = sim_name
        self.mix_time = mix_time
        self.output_dir = output_dir  # os.path.join(output_dir,sim_name)
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
        self.density = density
        self.activation_energy = activation_energy

        # below are default developer arguments which can be set through kwargs in sub classes for testing.
        self.legacy_bonding = False
        self.exclude_mixing_in_output = False # PLEASE NOTE THAT THE TRAJECTORY CHANGES WHEN THIS IS CHANGED!!
        self.init_file_name = os.path.join(self.output_dir, 'initial.hoomdxml')
        self.mixed_file_name = os.path.join(self.output_dir, 'mixed.hoomdxml')
        self.shrink_time = 1.0
        self.shrink = True
        self.ext_init_struct_path = None
        self.log_curing = False
        self.curing_log_period = 1e5
        self.curing_log = []
        self.bond_rank_log = []

        # for tests which compare simulation result against a benchmark
        # please see issue 6 for more details
        # (https://bitbucket.org/cmelab/getting-started/issues/6/different-trajectory-obtained-when-using).
        self.reset_random_after_initialize = False

    def get_sim_name(self):
        return self.simulation_name

    @abstractmethod
    def set_initial_structure(self):
        pass

    @abstractmethod
    def initialize_system_from_file(self, init_file_path=None, use_time_step_from_file=True):
        pass

    @abstractmethod
    def setup_mixing_run(self):
        pass

    @abstractmethod
    def setup_md_run(self):
        pass

    @staticmethod
    def initialize_context():
        try:
            __IPYTHON__
            run_from_ipython = True
        except NameError:
            run_from_ipython = False
        if run_from_ipython:
            print('Initializing HOOMD in ipython')
            hoomd.context.initialize('--mode=cpu')
        else:
            hoomd.context.initialize()

    @abstractmethod
    def get_curing_percentage(self, step):
        pass

    @abstractmethod
    def calculate_curing_percentage(self, step):
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
                 overwrite=True, static=['attribute'])
        deprecated.analyze.msd(groups=self.msd_groups, period=self.log_write, overwrite=True,
                               filename=os.path.join(self.output_dir, 'msd.log'), header_prefix='#')

    def run_mixing(self):
        md.integrate.mode_standard(dt=self.dt)
        md.integrate.nve(group=hoomd.group.all())
        hoomd.run(self.mix_time)
        if self.mixed_file_name.endswith('.hoomdxml'):
            deprecated.dump.xml(group=hoomd.group.all(), filename=self.mixed_file_name, all=True)
        elif self.mixed_file_name.endswith('.gsd'):
            hoomd.dump.gsd(group=hoomd.group.all(), filename=self.mixed_file_name, overwrite=True, period=None)

    def run_md(self):
        md.integrate.mode_standard(dt=self.dt)
        md.integrate.nve(group=hoomd.group.all())
        hoomd.run(self.md_time)
        deprecated.dump.xml(group=hoomd.group.all(),
                            filename=os.path.join(self.output_dir,
                                                  'final.hoomdxml'), all=True)

    def initialize(self):
        print('Initializing {}'.format(self.simulation_name))
        if not os.path.exists(self.output_dir):
            print('Creating simulation folder: {}'.format(self.output_dir))
            os.makedirs(self.output_dir)

        self.initialize_context()
        if self.ext_init_struct_path is None:
            self.set_initial_structure()
        else:
            print('Loading external initial structure : ', self.ext_init_struct_path)
            self.initialize_system_from_file(self.ext_init_struct_path)
            deprecated.dump.xml(group=hoomd.group.all(), filename=self.init_file_name, all=True)


        del self.system  # needed for re initializing hoomd
        self.initialize_context()
        self.initialize_system_from_file(self.init_file_name)
        deprecated.dump.xml(group=hoomd.group.all(), filename=self.init_file_name, all=True)

        self.setup_mixing_run()
        self.configure_outputs()
        self.run_mixing()

    def run(self):
        del self.system
        self.initialize_context()
        if self.exclude_mixing_in_output:
            use_time_step_from_file = False
        else:
            use_time_step_from_file = True
        self.initialize_system_from_file(self.mixed_file_name, use_time_step_from_file=use_time_step_from_file)

        print('Running MD for {}'.format(self.simulation_name))
        self.setup_md_run()
        self.configure_outputs()
        if self.bond is True:
            log = hoomd.analyze.log(filename=None, quantities=["temperature"], period=self.bond_period)
            if self.legacy_bonding is True:
                bonding_callback = LegacyBonding(system=self.system, groups=self.msd_groups, log=log,
                                                 activation_energy=self.activation_energy)
            else:
                bonding_callback = FreudBonding(system=self.system, groups=self.msd_groups, log=log,
                                                activation_energy=self.activation_energy)
            bond_callback = hoomd.analyze.callback(callback=bonding_callback, period=self.bond_period)

            if self.log_curing is True:
                curing_callback = hoomd.analyze.callback(callback=self.calculate_curing_percentage,
                                                         period=self.curing_log_period, phase=-1)
        hoomd.util.quiet_status()
        self.run_md()

        if self.bond is True:
            self.bond_rank_log = bonding_callback.bond_rank_log
            bond_callback.disable()
            log.disable()
            if self.log_curing is True:
                curing_callback.disable()

    def execute(self):
        print('Executing {}'.format(self.simulation_name))
        self.initialize()
        if self.reset_random_after_initialize:
            random.seed(12345)
        self.run()
        print("Finished executing {}".format(self.simulation_name))
