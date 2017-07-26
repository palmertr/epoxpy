import os
from epoxpy.simulation import Simulation
import hoomd
from hoomd import md
from hoomd import deprecated
from hoomd import dump
from abc import ABCMeta, abstractmethod
import numpy.random as rd
import numpy as np
import os, errno


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
                 bond=False, bond_period=1e1, box=[3, 3, 3], dt=1e-2, density=1.0, activation_energy=0.1,
                 sec_bond_weight=500.0, stop_bonding_after=None):
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
        self.sec_bond_weight = sec_bond_weight
        self.bonding = None
        self.bond_rank_hist_file = 'bond_rank_hist.log'
        self.log_bond_temp = None
        self.bond_callback = None
        self.dybond_updater = None
        self.stop_bonding_after = stop_bonding_after
        self.stop_dybond_updater_callback = None
        self.nl = None

        # below are default developer arguments which can be set through kwargs in sub classes for testing.
        self.nl_tuning = False
        self.profile_run = False
        self.legacy_bonding = False
        self.use_dybond_plugin = True
        self.exclude_mixing_in_output = False # PLEASE NOTE THAT THE TRAJECTORY CHANGES WHEN THIS IS CHANGED!!
        self.resume_file_name = os.path.join(self.output_dir, 'final.hoomdxml')
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

    @staticmethod
    def silent_remove(filename):
        try:
            os.remove(filename)
        except OSError as e:  # this would be "except OSError, e:" before Python 2.6
            if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
                raise  # re-raise exception if a different error occurred

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
        quantities=["pair_dpd_energy", "volume", "momentum", "potential_energy", "kinetic_energy", 
                    "temperature", "pressure", "bond_harmonic_energy"]
        if self.dybond_updater is not None:
            quantities.append("bond_percent(A-B)")
            #quantities.append("avg_num_failed_bonds")
            #quantities.append("avg_num_neighbors")

        print('quantitites being logged:',quantities)
        hoomd.analyze.log(filename=os.path.join(self.output_dir, 'out.log'),
                          quantities=quantities, period=self.log_write,
                          header_prefix='#', overwrite=True)
        dump.dcd(filename=os.path.join(self.output_dir, 'traj.dcd'), period=self.dcd_write, overwrite=True)
        dump.gsd(filename=os.path.join(self.output_dir, 'data.gsd'), period=self.dcd_write, group=hoomd.group.all(),
                 overwrite=True, static=['attribute'])
        deprecated.analyze.msd(groups=self.msd_groups, period=self.log_write, overwrite=True,
                               filename=os.path.join(self.output_dir, 'msd.log'), header_prefix='#')
        self.silent_remove(os.path.join(self.output_dir, self.bond_rank_hist_file))

    def run_mixing(self):
        md.integrate.mode_standard(dt=self.dt)
        md.integrate.nve(group=hoomd.group.all())
        hoomd.run(self.mix_time)
        if self.mixed_file_name.endswith('.hoomdxml'):
            deprecated.dump.xml(group=hoomd.group.all(), filename=self.mixed_file_name, all=True)
        elif self.mixed_file_name.endswith('.gsd'):
            hoomd.dump.gsd(group=hoomd.group.all(), filename=self.mixed_file_name, overwrite=True, period=None)

    @staticmethod
    def init_velocity(n, temp):
        v = rd.random((n, 3))
        v -= 0.5
        meanv = np.mean(v, 0)
        meanv2 = np.mean(v ** 2, 0)
        # fs = np.sqrt(3.*temp/meanv2)
        fs = np.sqrt(temp / meanv2)
        # print('scaling factor:{}'.format(fs))
        # print('v0:{}'.format(v))
        v = (v - meanv)  # shifts the average velocity of the simulation to 0
        v *= fs  # scaling velocity to match the desired temperature
        return v

    def initialize_snapshot_temperature(self, snapshot, temp):
        v = self.init_velocity(snapshot.particles.N, temp)
        snapshot.particles.velocity[:] = v[:]
        return snapshot

    def set_initial_particle_velocities(self, kT):
        snapshot = self.system.take_snapshot()
        snapshot = self.initialize_snapshot_temperature(snapshot, kT)
        self.system.restore_snapshot(snapshot)
        print('Reset the system temperature to {} kT after mixing'.format(kT))

    def run_md(self):
        md.integrate.mode_standard(dt=self.dt)
        md.integrate.nve(group=hoomd.group.all())
        first_target_temperature = self.temp_prof.temperature_profile[0][1]
        self.set_initial_particle_velocities(first_target_temperature)
        hoomd.run(self.md_time)
        if self.profile_run:
            hoomd.run(int(self.md_time*0.1),profile=True)# run 10% of the simulation time to calculate performance
        if self.nl_tuning:
            print('-----------------Disabling bonding and starting neighbourlist tuning-------------------')
            self.get_curing_percentage()
            if self.bond:
                if self.use_dybond_plugin:
                    self.dybond_updater.disable()
                else:
                    self.bond_callback.disable
            self.nl.tune(warmup=20000,
                         r_min=0.01,
                         r_max=2.00,
                         jumps=10,
                         steps=5000,
                         set_max_check_period=False)

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

        if self.log_curing is True and self.use_dybond_plugin is False:
            curing_callback = hoomd.analyze.callback(callback=self.calculate_curing_percentage,
                                                     period=self.curing_log_period, phase=-1)
        hoomd.util.quiet_status()
        self.run_md()

        if self.bond is True:
            if self.use_dybond_plugin is False:
                self.bond_callback.disable()
                if self.log_curing is True and self.use_dybond_plugin is False:
                    curing_callback.disable()
                self.log_bond_temp.disable()
    def execute(self):
        print('Executing {}'.format(self.simulation_name))
        self.initialize()
        if self.reset_random_after_initialize:
            random.seed(12345)
        self.run()
        print("Finished executing {}".format(self.simulation_name))
