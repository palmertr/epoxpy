from epoxy_simulation import EpoxySimulation
from epoxpy.lib import A, B, C, Epoxy_A_10_B_20_C10_2_Blend
import hoomd
from hoomd import md
from hoomd import deprecated
import mbuild as mb
import os


class ABCTypeEpoxySimulation(EpoxySimulation):
    """Simulations class for setting initial condition and force field specific to the ABC coarse grained Epoxy blend.
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
    def __init__(self, sim_name, mix_time, mix_kt, temp_prof, log_write=100, dcd_write=100, num_a=10, num_b=20,
                 num_c=2, n_mul=1.0, output_dir=os.getcwd(), bond=False,
                 bond_period=1e1, box=[3, 3, 3], dt=1e-2, **kwargs):
        EpoxySimulation.__init__(self, sim_name, mix_time=mix_time, mix_kt=mix_kt, temp_prof=temp_prof,
                                 log_write=log_write, dcd_write=dcd_write, output_dir=output_dir, bond=bond,
                                 bond_period=bond_period, box=box, dt=dt)
        self.num_a = num_a
        self.num_b = num_b
        self.num_c = num_c
        self.n_mul = n_mul
        self.dpd = None
        self.harmonic = None
        self.group_a = None
        self.group_b = None
        self.group_c = None

        print('kwargs passed into ABCTypeEpoxySimulation: {}'.format(kwargs))
        # setting developer variables through kwargs for testing.
        for key, value in kwargs.items():
            setattr(self, key, value)

    def set_initial_structure(self):
        blend = Epoxy_A_10_B_20_C10_2_Blend()
        mix_box = mb.packing.fill_box(blend, self.n_mul, box=self.box, overlap=0.050)
        file_name = os.path.join(self.output_dir, 'initial.gsd')
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

        if self.exclude_mixing_in_output is True:
            if self.init_file_name.endswith('.hoomdxml'):
                deprecated.dump.xml(group=hoomd.group.all(), filename=self.init_file_name, all=True)
            elif self.init_file_name.endswith('.gsd'):
                hoomd.dump.gsd(group=hoomd.group.all(), filename=self.init_file_name, overwrite=True, period=None)

            del self.system, snapshot  # needed for re initializing hoomd after randomize

    def setup_mixing_run(self):
        # Mix Step/MD Setup
        self.group_a = hoomd.group.type(name='a-particles', type='A')
        self.group_b = hoomd.group.type(name='b-particles', type='B')
        self.group_c = hoomd.group.type(name='c-particles', type='C')
        self.msd_groups = [self.group_a, self.group_b, self.group_c]

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

    def setup_md_run(self):
        self.group_a = hoomd.group.type(name='a-particles', type='A')
        self.group_b = hoomd.group.type(name='b-particles', type='B')
        self.group_c = hoomd.group.type(name='c-particles', type='C')
        self.msd_groups = [self.group_a, self.group_b, self.group_c]

        if self.exclude_mixing_in_output is True:
            nl = md.nlist.cell()
            profile = self.temp_prof.get_profile()
            print('temperature profile {}'.format(profile.points))
            #self.dpd.set_params(kT=profile)
            self.dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=profile, seed=0)
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
        else:
            profile = self.temp_prof.get_profile()
            print('temperature profile {}'.format(profile.points))
            self.dpd.set_params(kT=profile)
