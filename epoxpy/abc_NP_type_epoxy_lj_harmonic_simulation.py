from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
import hoomd
import cme_utils
from hoomd import md
import epoxpy.common as cmn
import mbuild as mb
import epoxpy.init as my_init
from epoxpy.lib import A, B, C, Sphere


class ABCNPTypeEpoxyLJHarmonicSimulation(ABCTypeEpoxySimulation):
    """Simulations class for ABCTypeEpoxySimulation where LJ is used as the
    conservative force and uses the langevin integrator.
       """
    def __init__(self,
                 sim_name,
                 mix_time,
                 mix_kt,
                 temp_prof,
                 AA_interaction=10.0,
                 AB_interaction=1.0,
                 AC_interaction=1.0,
                 BC_interaction=1.0,
                 AA_alpha=1.0,
                 AB_alpha=0.0,
                 AC_alpha=0.0,
                 BC_alpha=0.0,
                 num_spheres=10,
                 tau=0.1,
                 tauP=0.2,
                 P=1.0,
                 integrator=cmn.Integrators.LANGEVIN.name,
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
        self.AA_alpha = AA_alpha
        self.AB_alpha = AB_alpha
        self.AC_alpha = AC_alpha
        self.BC_alpha = BC_alpha
        self.num_spheres = num_spheres
        self.integrator = cmn.Integrators[integrator]
        self.tau = tau
        self.tauP = tauP
        self.P = P
        self._exclude_bonds_from_nlist = True

    def exclude_bonds_from_nlist(self):
        return self._exclude_bonds_from_nlist

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["pair_lj_energy", "bond_harmonic_energy"]
        return log_quantities

    def get_non_bonded_neighbourlist(self):
        nl = md.nlist.tree()  # cell()
        nl.reset_exclusions(exclusions=['bond']);
        return nl

    def set_initial_structure(self):
        print('========INITIAIZING STRUCTURE==========')
        num_beads_in_sphere = 65
        Sphere.mass = num_beads_in_sphere
        Sphere.name = 'Sphere'
        desired_box_volume = ((A.mass*self.num_a) + (B.mass*self.num_b) + (Sphere.mass*self.num_spheres)) / self.density
        desired_box_dim = (desired_box_volume ** (1./3.))
        reduced_density = self.density/10
        ex_box_vol = ((A.mass * self.num_a) + (B.mass * self.num_b) + (Sphere.mass*self.num_spheres)) / reduced_density
        expanded_box_dim = (ex_box_vol ** (1. / 3.))
        half_L = expanded_box_dim/2
        box = mb.Box(mins=[-half_L, -half_L, -half_L], maxs=[half_L, half_L, half_L])
        if self.old_init:
            raise NotImplementedError('Spherical Nano Particles not implemented in old init method')
        else:
            print("\n\n ===USING MBUILD INIT=== \n\n")
            if self.shrink is False:
                print('shrink=False is deprecated.')
            print('Packing {} A particles, {} B particles and {} C65s ..'.format(self.num_a,
                                                                                 self.num_b,
                                                                                 self.num_spheres))
            
            
            mix_box = mb.packing.fill_box([A(), B(), Sphere(n=num_beads_in_sphere)],
                                          [self.num_a, self.num_b, self.num_spheres],
                                          box=box)  # ,overlap=0.5)
            if self.num_spheres > 0: 
                mix_box.label_rigid_bodies(discrete_bodies = 'Sphere')

            if self.init_file_name.endswith('.hoomdxml'):
                mix_box.save(self.init_file_name, overwrite=True)
            elif self.init_file_name.endswith('.gsd'):
                mix_box.save(self.init_file_name, write_ff=False, overwrite=True)

            if self.init_file_name.endswith('.hoomdxml'):
                self.system = cme_utils.manip.convert_rigid.init_wrapper(xmlfile=self.init_file_name)
            elif self.init_file_name.endswith('.gsd'):
                self.system = hoomd.init.read_gsd(self.init_file_name)

            print('Initial box dimension: {}'.format(self.system.box.dimensions))


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
        #if self.num_b > 0 and self.num_c10 > 0:
        #    harmonic = md.bond.harmonic()
        #    harmonic.bond_coeff.set('C-C', k=self.CC_bond_const, r0=self.CC_bond_dist)
        #    harmonic.bond_coeff.set('A-B', k=self.AB_bond_const, r0=self.AB_bond_dist)
        lj = md.pair.lj(r_cut=2.5, nlist=self.nl)
        lj.pair_coeff.set('A', 'A', epsilon=self.AA_interaction, sigma=1.0, alpha=self.AA_alpha)
        lj.pair_coeff.set('B', 'B', epsilon=self.AA_interaction, sigma=1.0, alpha=self.AA_alpha)
        lj.pair_coeff.set('np', 'np', epsilon=self.AA_interaction, sigma=1.0, alpha=self.AA_alpha)

        lj.pair_coeff.set('A', 'B', epsilon=self.AB_interaction, sigma=1.0, alpha=self.AB_alpha)
        lj.pair_coeff.set('A', 'np', epsilon=self.AC_interaction, sigma=1.0, alpha=self.AC_alpha)
        lj.pair_coeff.set('B', 'np', epsilon=self.BC_interaction, sigma=1.0, alpha=self.BC_alpha)

    def setup_integrator(self, stage):
        print('=============Setting up {} integrator for {}'.format(self.integrator.name, stage.name))
        if stage == cmn.Stages.MIXING:
            temperature = self.mix_kT
            dt = self.mix_dt
            print('========= MIXING TEMPERATURE:', temperature, '=============')
        elif stage == cmn.Stages.CURING:
            profile = self.temp_prof.get_profile()
            temperature = profile
            dt = self.md_dt
            print('========= CURING TEMPERATURE:', temperature, '=============')
        md.integrate.mode_standard(dt=dt)
        rigid = hoomd.group.rigid_center()
        if self.integrator == cmn.Integrators.LANGEVIN:
            integrator = md.integrate.langevin(group=rigid,
                                               kT=temperature,
                                               seed=1223445,
                                               noiseless_t=False,
                                               noiseless_r=False)
            integrator.set_gamma('A', gamma=self.gamma)
            integrator.set_gamma('B', gamma=self.gamma)
            integrator.set_gamma('C', gamma=self.gamma)
        elif self.integrator == cmn.Integrators.NPT:
            integrator = md.integrate.npt(group=hoomd.group.all(),
                                          tau=self.tau,
                                          tauP=self.tauP,
                                          P=self.P,
                                          kT=temperature)
