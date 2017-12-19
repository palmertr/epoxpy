from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
from epoxpy.lib import A, B, C, C10, Epoxy_A_10_B_20_C10_2_Blend
import epoxpy.init as my_init
import mbuild as mb
from hoomd import deprecated
import hoomd
from hoomd import md
import epoxpy.common as cmn


class ABCTypeEpoxyDPDSimulation(ABCTypeEpoxySimulation):
    """Simulations class for ABCTypeEpoxySimulation where DPD is used as the
    conservative force.
       """
    def __init__(self,
                 sim_name,
                 mix_time,
                 mix_kt,
                 temp_prof,
                 AA_interaction=25.0,
                 AB_interaction=35.0,
                 AC_interaction=35.0,
                 BC_interaction=35.0,
                 *args,
                 **kwargs):
        ABCTypeEpoxySimulation.__init__(self,
                                        sim_name,
                                        mix_time,
                                        mix_kt,
                                        temp_prof,
                                        *args,
                                        **kwargs)
        self.AA_interaction = AA_interaction
        self.AB_interaction = AB_interaction
        self.AC_interaction = AC_interaction
        self.BC_interaction = BC_interaction
        self.dpd = None

    def set_initial_structure(self):
        print('========INITIAIZING FOR DPD==========')
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

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["pair_dpd_energy","bond_harmonic_energy"]
        return log_quantities

    def get_non_bonded_neighbourlist(self):
        nl = md.nlist.tree()  # cell()
        nl.reset_exclusions(exclusions=[]);
        return nl

    def setup_force_fields(self, stage):
        if self.DEBUG:
            print('=============force fields parameters==============')
            print('self.CC_bond_const', self.CC_bond_const)
            print('self.CC_bond_dist', self.CC_bond_dist)
            print('self.AB_bond_const', self.AB_bond_const)
            print('self.AB_bond_dist', self.AB_bond_dist)
            print('self.AA_interaction', self.AA_interaction)
            print('self.AB_interaction', self.AB_interaction)
            print('self.AC_interaction', self.AC_interaction)
            print('self.BC_interaction', self.BC_interaction)
            print('self.gamma', self.gamma)

        if stage == cmn.Stages.MIXING:
            temperature = self.mix_kT
            print('========= MIXING TEMPERATURE:', temperature, '=============')
        elif stage == cmn.Stages.CURING:
            profile = self.temp_prof.get_profile()
            temperature = profile
            print('========= CURING TEMPERATURE:', temperature, '=============')

        if self.num_b > 0 and self.num_c10 > 0:
            harmonic = md.bond.harmonic()
            harmonic.bond_coeff.set('C-C', k=self.CC_bond_const, r0=self.CC_bond_dist)
            harmonic.bond_coeff.set('A-B', k=self.AB_bond_const, r0=self.AB_bond_dist)
        dpd = md.pair.dpd(r_cut=1.0, nlist=self.nl, kT=temperature, seed=123456)
        dpd.pair_coeff.set('A', 'A', A=self.AA_interaction, gamma=self.gamma)
        dpd.pair_coeff.set('B', 'B', A=self.AA_interaction, gamma=self.gamma)
        dpd.pair_coeff.set('C', 'C', A=self.AA_interaction, gamma=self.gamma)

        dpd.pair_coeff.set('A', 'B', A=self.AB_interaction, gamma=self.gamma)
        dpd.pair_coeff.set('A', 'C', A=self.AC_interaction, gamma=self.gamma)
        dpd.pair_coeff.set('B', 'C', A=self.BC_interaction, gamma=self.gamma)

    def setup_integrator(self, stage):
        if stage == cmn.Stages.MIXING:
            dt = self.mix_dt
            print('========= MIXING dt:', dt, '=============')
        elif stage == cmn.Stages.CURING:
            dt = self.md_dt
            print('========= CURING dt:', dt, '=============')
        md.integrate.mode_standard(dt=dt)
        md.integrate.nve(group=hoomd.group.all())
