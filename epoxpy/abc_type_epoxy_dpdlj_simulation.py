from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
from epoxpy.lib import A, B, C, C10, Epoxy_A_10_B_20_C10_2_Blend
import epoxpy.init as my_init
import hoomd
import hoomd.dybond_plugin as db
from hoomd import md
from hoomd import deprecated
from hoomd import variant
import mbuild as mb

class ABCTypeEpoxyDPDLJSimulation(ABCTypeEpoxySimulation):
    """Simulations class for ABCTypeEpoxySimulation where LJ is used as the
    conservative force in the DPD force field.
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
                 shrink_time=1e6,
                 AA_alpha=1.0,
                 AB_alpha=0.0,
                 AC_alpha=0.0,
                 BC_alpha=0.0,
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
        self.shrink_time = shrink_time
        self.shrinkT = 5.0

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["pair_dpdlj_energy"]
        return log_quantities

    def set_initial_structure(self):
        print('========INITIAIZING FOR LJ==========')
        desired_box_volume = ((A.mass*self.num_a) + (B.mass*self.num_b) + (C10.mass*self.num_c10)) / self.density
        desired_box_dim = (desired_box_volume ** (1./3.))
        self.box = [desired_box_dim, desired_box_dim,
                    desired_box_dim]
        if self.old_init == True:
            print("\n\n ===USING OLD INIT=== \n\n")
            As = my_init.Bead(btype="A", mass=A.mass)
            Bs = my_init.Bead(btype="B", mass=B.mass)
            C10s = my_init.PolyBead(btype="C", mass = 1.0, N = 10) # Hardcode C10, with mon-mass 1.0
            snap = my_init.init_system({As : int(self.num_a), Bs
                                        :int(self.num_b), C10s :
                                        int(self.num_c10)}, self.density/10)
            self.system = hoomd.init.read_snapshot(snap)
        else:
            if self.shrink is True:
                print('Packing {} A particles ..'.format(self.num_a))
                mix_box = mb.packing.fill_box(A(), self.num_a,
                                              box=self.box)#,overlap=0.5)
                mix_box = mb.packing.solvate(mix_box, B(), self.num_b,
                                             box=self.box)#,overlap=0.5)
                print('Packing {} B particles ..'.format(self.num_b))
                mix_box = mb.packing.solvate(mix_box, C10(), self.num_c10,
                                             box=self.box)#,overlap=0.5)
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
            self.setup_mixing_run()
            size_variant =\
            variant.linear_interp([(0,self.system.box.Lx),(self.shrink_time,desired_box_dim)])
            md.integrate.mode_standard(dt=self.dt)
            md.integrate.langevin(group=hoomd.group.all(),
                                  kT=self.shrinkT,
                                  seed=1223445)#self.seed)
            resize=hoomd.update.box_resize(L=size_variant)
            #hoomd.update.box_resize(period=1, L=desired_box_dim)
            hoomd.run(self.shrink_time)
            snapshot = self.system.take_snapshot()
            print('Initial box dimension: {}'.format(snapshot.box))

        if self.init_file_name.endswith('.hoomdxml'):
            deprecated.dump.xml(group=hoomd.group.all(), filename=self.init_file_name, all=True)
        elif self.init_file_name.endswith('.gsd'):
            hoomd.dump.gsd(group=hoomd.group.all(), filename=self.init_file_name, overwrite=True, period=None)

    def setup_mixing_run(self):
        # Mix Step/MD Setup
        super().setup_mixing_run()
        dpdlj = md.pair.dpdlj(r_cut=2.5, nlist=self.nl, kT=self.mix_kT, seed=123456)
        dpdlj.pair_coeff.set('A', 'A', epsilon=self.AA_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AA_alpha)
        dpdlj.pair_coeff.set('B', 'B', epsilon=self.AA_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AA_alpha)
        dpdlj.pair_coeff.set('C', 'C', epsilon=self.AA_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AA_alpha)

        dpdlj.pair_coeff.set('A', 'B', epsilon=self.AB_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AB_alpha)
        dpdlj.pair_coeff.set('A', 'C', epsilon=self.AC_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AC_alpha)
        dpdlj.pair_coeff.set('B', 'C', epsilon=self.BC_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.BC_alpha)

    def setup_md_run(self):
        super().setup_md_run()
        profile = self.temp_prof.get_profile()
        print('temperature profile {}'.format(profile.points))

        dpdlj = md.pair.dpdlj(r_cut=2.5, nlist=self.nl, kT=profile, seed=123456)
        dpdlj.set_params(kT=profile)
        dpdlj.pair_coeff.set('A', 'A', epsilon=self.AA_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AA_alpha)
        dpdlj.pair_coeff.set('B', 'B', epsilon=self.AA_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AA_alpha)
        dpdlj.pair_coeff.set('C', 'C', epsilon=self.AA_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AA_alpha)
                                                                                                          
        dpdlj.pair_coeff.set('A', 'B', epsilon=self.AB_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AB_alpha)
        dpdlj.pair_coeff.set('A', 'C', epsilon=self.AC_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.AC_alpha)
        dpdlj.pair_coeff.set('B', 'C', epsilon=self.BC_interaction, sigma=1.0 , gamma=self.gamma,alpha=self.BC_alpha)

