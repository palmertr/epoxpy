from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
from epoxpy.lib import A, B, C, C10, Epoxy_A_10_B_20_C10_2_Blend
import epoxpy.init as my_init
import hoomd
import hoomd.dybond_plugin as db
from hoomd import md
from hoomd import deprecated
from hoomd import variant
import mbuild as mb

class ABCTypeEpoxyNPTQuenchSimulation(ABCTypeEpoxySimulation):
    """Simulations class for ABCTypeEpoxySimulation where LJ is used as the
    conservative force and uses the npt integrator.
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
                 tau=0.1,
                 tauP=0.2,
                 P=1.0,
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
        self.tau = tau
        self.tauP = tauP
        self.P = P

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["bond_harmonic_energy"]
        return log_quantities

    def set_initial_structure(self):
        print('=================== WARNING!! set_initial_structure function\
              should not be called for quench simulations')

    def setup_forcefields(self):
        if self.DEBUG:
            print('=============force fields parameters==============')
            print('self.CC_bond_const',self.CC_bond_const)
            print('self.CC_bond_dist',self.CC_bond_dist)
            print('self.AB_bond_const',self.AB_bond_const)
            print('self.AB_bond_dist',self.AB_bond_dist)
            print('self.AA_interaction',self.AA_interaction)
            print('self.AB_interaction',self.AB_interaction)
            print('self.AC_interaction',self.AC_interaction)
            print('self.BC_interaction',self.BC_interaction)
            print('self.gamma',self.gamma)
        if self.num_b > 0 and self.num_c10 > 0:
            harmonic = md.bond.harmonic()
            harmonic.bond_coeff.set('C-C', k=self.CC_bond_const, r0=self.CC_bond_dist)
            harmonic.bond_coeff.set('A-B', k=self.AB_bond_const, r0=self.AB_bond_dist)
        lj = md.pair.lj(r_cut=2.5, nlist=self.nl)
        lj.pair_coeff.set('A', 'A', epsilon=self.AA_interaction, sigma=1.0 , alpha=self.AA_alpha)
        lj.pair_coeff.set('B', 'B', epsilon=self.AA_interaction, sigma=1.0 , alpha=self.AA_alpha)
        lj.pair_coeff.set('C', 'C', epsilon=self.AA_interaction, sigma=1.0 , alpha=self.AA_alpha)

        lj.pair_coeff.set('A', 'B', epsilon=self.AB_interaction, sigma=1.0 , alpha=self.AB_alpha)
        lj.pair_coeff.set('A', 'C', epsilon=self.AC_interaction, sigma=1.0 , alpha=self.AC_alpha)
        lj.pair_coeff.set('B', 'C', epsilon=self.BC_interaction, sigma=1.0 , alpha=self.BC_alpha)

    def setup_mixing_run(self):
        # Mix Step/MD Setup
        super().setup_mixing_run()
        self.nl.reset_exclusions(exclusions = ['bond']);
        self.setup_forcefields()
        bd=md.integrate.langevin(group=hoomd.group.all(),
                        kT=self.mix_kT,
                        seed=1223445,noiseless_t=False, noiseless_r=False)
        bd.set_gamma('A', gamma=self.gamma)
        bd.set_gamma('B', gamma=self.gamma)
        bd.set_gamma('C', gamma=self.gamma)

    def setup_md_run(self):
        super().setup_md_run()
        self.nl.reset_exclusions(exclusions = ['bond']);
        profile = self.temp_prof.get_profile()
        print('temperature profile {}'.format(profile.points))
        self.setup_forcefields()
        print('===========================md_dt=====================',self.md_dt)
        md.integrate.mode_standard(dt=self.md_dt)
        npt=md.integrate.npt(group=hoomd.group.all(),
                            tau=self.tau,
                            tauP=self.tauP,
                            P=self.P,
                            kT=profile)
