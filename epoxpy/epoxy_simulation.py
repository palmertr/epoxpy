import os
from simulation import Simulation
from bonding import Bonding
import hoomd
from hoomd import md
from hoomd import deprecated
from hoomd import dump
from epoxpy.lib import A, B, C, Epoxy_A_10_B_20_C10_2_Blend
import mbuild as mb


class EpoxySimulation(Simulation):
    """Simulations class for setting initial condition, force field and choosing MD engine.
       This simulation consists of three particle types (A, B and C). A, B and C particles are created in the
       ratio 10, 20 and 2 by default
       The force fields used are Dissipative Particle Dynamics (DPD) and Harmonic potential

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
       num_a    : number of A particles created.
       num_b    : number of B particles created.
       num_C    : number of C particles created.
       n_mul    : multiplying factor for the number of A, B and C particles to be created.
       output_dir: default is the working directory
       bond     : boolean value denoting whether to run the bonding routine for A's and B's
       bond_period: time interval between calls to the bonding routine
    """
    engine_name = 'HOOMD'

    def __init__(self, sim_name, mix_time, mix_kt, temp_prof, log_write=100, dcd_write=100, num_a=10, num_b=20,
                 num_c=2, n_mul=1.0, output_dir=os.getcwd(), bond=False, bond_period=1e1):
        Simulation.__init__(self, self.engine_name)
        self.simulation_name = sim_name
        self.n_mult = n_mul
        self.mix_time = mix_time
        self.output_dir = output_dir+"/"+sim_name+"/"
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
        self.num_a = num_a
        self.num_b = num_b
        self.num_c = num_c
        self.dpd = None
        self.harmonic = None
        self.group_a = None
        self.group_b = None
        self.group_c = None
        self.system = None

    def get_sim_name(self):
        return self.simulation_name

    def set_initial_structure(self):
        blend = Epoxy_A_10_B_20_C10_2_Blend()
        mix_box = mb.packing.fill_box(blend, self.n_mult, box=[3, 3, 3])
        file_name = 'initial.gsd'
        if file_name.endswith('.hoomdxml'):
            mix_box.save(file_name)
        elif file_name.endswith('.gsd'):
            mix_box.save(file_name, write_ff=False)

        if file_name.endswith('.hoomdxml'):
            self.system = hoomd.deprecated.init.read_xml(file_name)
        elif file_name.endswith('.gsd'):
            self.system = hoomd.init.read_gsd(file_name)

        print('Initial box dimension: {}'.format(self.system.box.dimensions))

        snapshot = self.system.take_snapshot(bonds=True)
        for p_id in range(snapshot.particles.N):
            p_types = snapshot.particles.types
            p_type = p_types[snapshot.particles.typeid[p_id]]
            if p_type == 'A':
                snapshot.particles.mass[p_id] = A.mass
            if p_type == 'B':
                snapshot.particles.mass[p_id] = B.mass
            if p_type == 'C':
                snapshot.particles.mass[p_id] = C.mass
        print(snapshot.bonds.types)
        snapshot.bonds.types = ['C-C', 'A-B']
        self.system.restore_snapshot(snapshot)
        self.system.bonds.add('A-B', 1, 15)

    def configure_outputs(self):
        print('Configuring outputs. output_dir: {}'.format(self.output_dir))
        print('log_write: {} dcd_write: {}'.format(self.log_write, self.dcd_write))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        hoomd.meta.dump_metadata(filename=self.output_dir + "metadata.json", indent=2)
        deprecated.dump.xml(group=hoomd.group.all(), filename=self.output_dir + "start.hoomdxml", all=True)
        hoomd.analyze.log(filename=self.output_dir + "out.log",
                          quantities=["pair_dpd_energy", "volume", "momentum", "potential_energy", "kinetic_energy",
                                      "temperature", "pressure", "bond_harmonic_energy"], period=self.log_write,
                          header_prefix='#', overwrite=True)
        dump.dcd(filename=self.output_dir + "traj.dcd", period=self.dcd_write, overwrite=True)
        dump.gsd(filename=self.output_dir + "data.gsd", period=self.dcd_write, group=hoomd.group.all(), overwrite=True,
                 static=['attribute', 'topology'])

    def initialize(self):
        print('Initializing {}'.format(self.simulation_name))
        #self.engine.set_initial_structure(self.get_initial_structure())
        self.set_initial_structure()

        # Mix Step/MD Setup
        self.group_a = hoomd.group.type(name='a-particles', type='A')
        self.group_b = hoomd.group.type(name='b-particles', type='B')
        self.group_c = hoomd.group.type(name='c-particles', type='C')

        nl = md.nlist.cell()
        self.dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=self.mix_kT, seed=0)
        # dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=1.0, gamma = 1.0)
        self.dpd.pair_coeff.set('A', 'A', A=1.0, gamma=1.0)
        self.dpd.pair_coeff.set('B', 'B', A=1.0, gamma=1.0)
        self.dpd.pair_coeff.set('C', 'C', A=1.0, gamma=1.0)

        self.dpd.pair_coeff.set('A', 'B', A=10.0, gamma=1.0)
        self.dpd.pair_coeff.set('A', 'C', A=10.0, gamma=1.0)
        self.dpd.pair_coeff.set('B', 'C', A=10.0, gamma=1.0)

        self.harmonic = md.bond.harmonic()
        # 0 for C-C bond. We know its assigned a bond id 0 from the hoomdxml file saved by mbuild
        self.harmonic.bond_coeff.set('C-C', k=100.0, r0=1.0)
        self.harmonic.bond_coeff.set('A-B', k=100.0, r0=1.0)

        self.configure_outputs()
        self.engine.run_initial_structure(self.mix_time, self.output_dir)

    def run(self):
        print('Running {}'.format(self.simulation_name))
        print('temperature profile {}'.format(self.temp_prof.get_profile()))
        self.dpd.set_params(kT=self.temp_prof.get_profile())
        deprecated.analyze.msd(groups=[self.group_a, self.group_b, self.group_c], period=self.log_write,
                               filename=self.output_dir + "msd.log", header_prefix='#')

        if self.bond is True:
            log = hoomd.analyze.log(filename=None, quantities=["temperature"], period=self.bond_period)
            bonding_callback = Bonding(system=self.system, epoxy_sim=self, log=log)
            bond_callback = hoomd.analyze.callback(callback=bonding_callback, period=self.bond_period)

        self.engine.run_md(self.md_time, self.output_dir)

        if self.bond is True:
            bond_callback.disable()
            log.disable()

    def output(self):
        import numpy as np
        data = np.genfromtxt(fname=self.output_dir + "out.log",skip_header=True)

    def execute(self):
        print('Executing {}'.format(self.simulation_name))
        self.initialize()
        self.run()
        self.output()
        print("sim fin")
