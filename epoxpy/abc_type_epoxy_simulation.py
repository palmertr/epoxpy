from epoxpy.epoxy_simulation import EpoxySimulation
from epoxpy.lib import A, B, C, C10, Epoxy_A_10_B_20_C10_2_Blend
from epoxpy.bonding import LegacyBonding, FreudBonding
import hoomd
import hoomd.dybond_plugin as db
from hoomd import md
from hoomd import deprecated
import mbuild as mb
import os
import numpy as np
import epoxpy.init as my_init
from collections import Counter


class ABCTypeEpoxySimulation(EpoxySimulation):
    MAX_A_BONDS=4
    MAX_B_BONDS=2#This ratio is stoichiometric by default (see num_a and num_b defaults)
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
    def __init__(self,
                 sim_name,
                 mix_time,
                 mix_kt,
                 temp_prof,
                 num_a=10,
                 num_b=20,
                 num_c10=2,
                 n_mul=1.0,
                 gamma=4.5,
                 stop_after_percent=100.0,
                 percent_bonds_per_step=0.0025,
                 AB_bond_const=100,
                 CC_bond_const=100,
                 AB_bond_dist=1,
                 CC_bond_dist=1,
                 *args,
                 **kwargs):
        EpoxySimulation.__init__(self,
                                 sim_name,
                                 mix_time,
                                 mix_kt,
                                 temp_prof,
                                 *args,
                                 **kwargs)
        self.num_a = int(num_a * n_mul)
        self.num_b = int(num_b * n_mul)
        self.num_c10 = int(num_c10 * n_mul)
        self.n_mul = n_mul
        self.harmonic = None
        self.group_a = None
        self.group_b = None
        self.group_c = None
        self.AB_bond_const=AB_bond_const
        self.CC_bond_const=CC_bond_const
        self.AB_bond_dist=AB_bond_dist
        self.CC_bond_dist=CC_bond_dist
        self.gamma = gamma
        self.stop_after_percent = stop_after_percent
        self.percent_bonds_per_step = percent_bonds_per_step
        print('kwargs passed into ABCTypeEpoxySimulation: {}'.format(kwargs))
        # setting developer variables through kwargs for testing.
        for key, value in kwargs.items():
            setattr(self, key, value)

    def set_initial_structure(self): # TODO this needs to be set in the base class or something
        desired_box_volume = ((A.mass*self.num_a) + (B.mass*self.num_b) + (C10.mass*self.num_c10)) / self.density
        desired_box_dim = (desired_box_volume ** (1./3.))
        self.box = [desired_box_dim, desired_box_dim, desired_box_dim]
        if self.old_init == True:
            print("\n\n ===USING OLD INIT=== \n\n")
            As = my_init.Bead(btype="A", mass=A.mass)
            Bs = my_init.Bead(btype="B", mass=B.mass)
            C10s = my_init.PolyBead(btype="C", mass = 1.0, N = 10) # Hardcode C10, with mon-mass 1.0
            snap = my_init.init_system({As : int(self.num_a), Bs :int(self.num_b), C10s : int(self.num_c10)}, self.density)
            self.system = hoomd.init.read_snapshot(snap)
        else:
            if self.shrink is True:
                print('Packing {} A particles ..'.format(self.num_a))
                mix_box = mb.packing.fill_box(A(), self.num_a, box=self.box)
                mix_box = mb.packing.solvate(mix_box, B(), self.num_b, box=self.box)
                print('Packing {} B particles ..'.format(self.num_b))
                mix_box = mb.packing.solvate(mix_box, C10(), self.num_c10, box=self.box)
                print('Packing {} C10 particles ..'.format(self.num_c10))
            else:
                blend = Epoxy_A_10_B_20_C10_2_Blend()
                mix_box = mb.packing.fill_box(blend, self.n_mul, box=self.box, overlap=0.050)

            if self.init_file_name.endswith('.hoomdxml'):
                mix_box.save(self.init_file_name)
            elif self.init_file_name.endswith('.gsd'):
                mix_box.save(self.init_file_name, write_ff=False)

            if self.init_file_name.endswith('.hoomdxml'):
                self.system = hoomd.deprecated.init.read_xml(self.init_file_name)
            elif self.init_file_name.endswith('.gsd'):
                self.system = hoomd.init.read_gsd(self.init_file_name)

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

        if self.shrink is True:
            hoomd.update.box_resize(period=1, L=desired_box_dim)
            hoomd.run(self.shrink_time)
            snapshot = self.system.take_snapshot()
            print('Initial box dimension: {}'.format(snapshot.box))

        if self.init_file_name.endswith('.hoomdxml'):
            deprecated.dump.xml(group=hoomd.group.all(), filename=self.init_file_name, all=True)
        elif self.init_file_name.endswith('.gsd'):
            hoomd.dump.gsd(group=hoomd.group.all(), filename=self.init_file_name, overwrite=True, period=None)

    def initialize_system_from_file(self, file_path, use_time_step_from_file=True):
        if use_time_step_from_file:
            time_step = None
        else:
            time_step = 0  # start simulation from start.

        if file_path.endswith('.hoomdxml'):
            self.system = hoomd.deprecated.init.read_xml(file_path, time_step=time_step)
        elif file_path.endswith('.gsd'):
            raise ValueError('Reading the most recent frame from gsd file is not yet implemented!')
            self.system = hoomd.init.read_gsd(file_path, frame=0, time_step=time_step)
        else:
            raise ValueError('No such file as {} exist on disk!'.format(file_path))
        snapshot = self.system.take_snapshot(bonds=True)
        snapshot.bonds.types = ['C-C', 'A-B']
        self.system.restore_snapshot(snapshot)

    def setup_mixing_run(self):
        # Mix Step/MD Setup
        self.group_a = hoomd.group.type(name='a-particles', type='A')
        self.group_b = hoomd.group.type(name='b-particles', type='B')
        self.group_c = hoomd.group.type(name='c-particles', type='C')
        self.msd_groups = [self.group_a, self.group_b, self.group_c]

        self.nl = md.nlist.cell()
        if self.num_b > 0 and self.num_c10 > 0:
            self.harmonic = md.bond.harmonic()
            self.harmonic.bond_coeff.set('C-C', k=self.CC_bond_const, r0=self.CC_bond_dist)
            self.harmonic.bond_coeff.set('A-B', k=self.AB_bond_const, r0=self.AB_bond_dist)
        self.nl.reset_exclusions(exclusions = []);

    def print_curing_and_stop_updater(self, bond_percent):
        print("HIT OUR TARGET: {}".format(self.stop_after_percent))
        print("PERCENT REACHED: {}".format(bond_percent))
        #self.dybond_updater.disable()# this causes segfault

    def stop_dybond_updater(self, timestep):
        if self.stop_dybond_updater_callback is not None:
            self.dybond_updater.disable() # first stop the updater
            self.stop_dybond_updater_callback.disable() # now stop the callback.
        else:
            hoomd.context.msg.warning('Call back for stopping the bonding is not set!')

    def setup_md_run(self):
        self.group_a = hoomd.group.type(name='a-particles', type='A')
        self.group_b = hoomd.group.type(name='b-particles', type='B')
        self.group_c = hoomd.group.type(name='c-particles', type='C')
        self.msd_groups = [self.group_a, self.group_b, self.group_c]

        self.nl = md.nlist.cell()

        if self.num_b > 0 and self.num_c10 > 0:
            self.harmonic = md.bond.harmonic()
            self.harmonic.bond_coeff.set('C-C', k=self.CC_bond_const, r0=self.CC_bond_dist)
            self.harmonic.bond_coeff.set('A-B', k=self.AB_bond_const, r0=self.AB_bond_dist)
            print('C-C bond constant:',self.CC_bond_const,'C-C bond dist:',self.CC_bond_dist)
        self.nl.reset_exclusions(exclusions = []);

        if self.bond is True:
            if self.use_dybond_plugin is True:
                self.dybond_updater = db.update.dybond(self.nl, group=hoomd.group.all(), period=self.bond_period)
                self.dybond_updater.set_params(bond_type='A-B',A='A',
                                               A_fun_groups=ABCTypeEpoxySimulation.MAX_A_BONDS,B='B',
                                               B_fun_groups=ABCTypeEpoxySimulation.MAX_B_BONDS,
                                               Ea=self.activation_energy,
                                               rcut=self.bond_radius,alpha=self.sec_bond_weight,
                                               percent_bonds_per_step=self.percent_bonds_per_step,
                                               stop_after_percent=self.stop_after_percent,
                                               callback=self.print_curing_and_stop_updater)

                if self.stop_bonding_after is not None:
                    self.stop_dybond_updater_callback = hoomd.analyze.callback(callback=self.stop_dybond_updater,
                                                                           period=self.stop_bonding_after)
            else:
                self.log_bond_temp = hoomd.analyze.log(filename=None, quantities=["temperature"],
                                                        period=self.bond_period)
                if self.legacy_bonding is True:
                    self.bonding = LegacyBonding(system=self.system, groups=self.msd_groups, log=self.log_bond_temp,
                                                 activation_energy=self.activation_energy,
                                             sec_bond_weight=self.sec_bond_weight)
                else:
                    self.bonding = FreudBonding(system=self.system, groups=self.msd_groups, log=self.log_bond_temp,
                                                activation_energy=self.activation_energy,
                                            sec_bond_weight=self.sec_bond_weight)
                self.bond_callback = hoomd.analyze.callback(callback=self.bonding, period=self.bond_period)

    def total_possible_bonds(self):
        if self.num_b * FreudBonding.MAX_B_BONDS > self.num_a * FreudBonding.MAX_A_BONDS:
           possible_bonds = (self.num_a * FreudBonding.MAX_A_BONDS)
        else:
            possible_bonds = (self.num_b * FreudBonding.MAX_B_BONDS)
        return possible_bonds

    def get_curing_percentage(self):
        n_bonds = 0
        if self.system is not None:
            snapshot = self.system.take_snapshot(bonds=True)
            n_bonds = len(snapshot.bonds.group) - (self.num_c10 * 9)
        possible_bonds = self.total_possible_bonds()
        bond_percent = (n_bonds / possible_bonds) * 100.
        #print('possible bonds:{}, bonds made:{}, cure percent: {}'.format(possible_bonds, n_bonds, bond_percent))
        return bond_percent

    '''deprecated: Used with freud and legacy bonding
    '''
    def calculate_curing_percentage(self, step):
        bond_percent = self.get_curing_percentage()
        self.curing_log.append((step, bond_percent))

        group_a_idx = []
        for p in self.group_a:
            group_a_idx.append(p.tag)
        #print(group_a_idx)
        group_b_idx = []
        for p in self.group_b:
            group_b_idx.append(p.tag)
        #print(group_b_idx)

        dic_vals = [self.bonding.rank_dict.get(k, 0) for k in group_a_idx]
        keys = (list(Counter(dic_vals).keys()))
        values = (list(Counter(dic_vals).values()))
        #print(keys)
        #print(values)
        row = np.zeros(5)
        for i in range(0, 4):
            if i + 1 in keys:
                row[i] = values[keys.index(i + 1)]
        row = [row]
        this_row = row[0]
        p_bonds = (100. * this_row[0]) / self.num_a
        s_bonds = (100. * this_row[1]) / self.num_a
        t_bonds = (100. * this_row[2]) / self.num_a
        q_bonds = (100. * this_row[3]) / self.num_a

        dic_vals = [self.bonding.rank_dict.get(k, 0) for k in group_b_idx]
        keys = (list(Counter(dic_vals).keys()))
        values = (list(Counter(dic_vals).values()))
        #print(keys)
        #print(values)
        if 0 in keys:
            primary_b = values[keys.index(0)]
        else:
            primary_b = 0
        primary_b = (100. * primary_b) / self.num_b

        row = [(step, bond_percent, p_bonds, s_bonds, t_bonds, q_bonds, primary_b)]
        with open(os.path.join(self.output_dir, self.bond_rank_hist_file), 'ab') as f_handle:
            np.savetxt(f_handle, row)
