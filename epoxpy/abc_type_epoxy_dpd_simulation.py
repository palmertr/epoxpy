from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
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
                 nlist=cmn.NeighbourList.CELL.name,
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
        self.nlist = cmn.NeighbourList[nlist]
        self.dpd = None
        self._exclude_bonds_from_nlist = False

    def exclude_bonds_from_nlist(self):
        return self._exclude_bonds_from_nlist

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["pair_dpd_energy","bond_harmonic_energy"]
        return log_quantities

    def get_non_bonded_neighbourlist(self):
        if self.nlist == cmn.NeighbourList.CELL:
            nl = md.nlist.cell()
        else:
            nl = md.nlist.tree()
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
