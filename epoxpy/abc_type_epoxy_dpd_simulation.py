from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
from hoomd import md

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

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["pair_dpd_energy"]
        return log_quantities

    def setup_mixing_run(self):
        # Mix Step/MD Setup
        super().setup_mixing_run()
        self.dpd = md.pair.dpd(r_cut=1.0, nlist=self.nl, kT=self.mix_kT, seed=123456)
        self.dpd.pair_coeff.set('A', 'A', A=self.AA_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('B', 'B', A=self.AA_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('C', 'C', A=self.AA_interaction, gamma=self.gamma)

        self.dpd.pair_coeff.set('A', 'B', A=self.AB_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('A', 'C', A=self.AC_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('B', 'C', A=self.BC_interaction, gamma=self.gamma)

    def setup_md_run(self):
        super().setup_md_run()
        profile = self.temp_prof.get_profile()
        print('temperature profile {}'.format(profile.points))
        self.dpd = md.pair.dpd(r_cut=1.0, nlist=self.nl, kT=profile, seed=123456)
        self.dpd.set_params(kT=profile)
        self.dpd.pair_coeff.set('A', 'A', A=self.AA_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('B', 'B', A=self.AA_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('C', 'C', A=self.AA_interaction, gamma=self.gamma)

        self.dpd.pair_coeff.set('A', 'B', A=self.AB_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('A', 'C', A=self.AC_interaction, gamma=self.gamma)
        self.dpd.pair_coeff.set('B', 'C', A=self.BC_interaction, gamma=self.gamma)

