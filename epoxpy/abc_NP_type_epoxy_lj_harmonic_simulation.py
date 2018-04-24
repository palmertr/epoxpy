from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation
import hoomd
import cme_utils
from hoomd import md
import epoxpy.common as cmn
import mbuild as mb
import os
from epoxpy.lib import A, B, AB, Sphere
from hoomd import deprecated
import cme_utils.manip.convert_rigid as init_rigid
from hoomd import variant
import numpy
from numpy import pi


class ABCNPTypeEpoxyLJHarmonicSimulation(ABCTypeEpoxySimulation):
    """Simulations class for ABCTypeEpoxySimulation where LJ is used as the
    conservative force and uses the langevin integrator.
       """
    def __init__(self,
                 sim_name,
                 mix_time,
                 mix_kt,
                 temp_prof,
                 AA_interaction=1.0,
                 AB_interaction=1.0,
                 AC_interaction=1.0,
                 BC_interaction=1.0,
                 Anp_interaction=1.0,
                 Bnp_interaction=1.0,
                 Arc_interaction=1.0,
                 Brc_interaction=1.0,
                 Crc_interaction=1.0,
                 npnp_interaction=1.0,
                 nprc_interaction=1.0,
                 rcrc_interaction=1.0,
                 Cnp_interaction=1.0,
                 AA_alpha=1.0,
                 AB_alpha=1.0,
                 AC_alpha=1.0,
                 BC_alpha=1.0,
                 num_spheres=5,
                 radius_sphere=5,
                 radius_a=0.5,
                 radius_b=0.5,
                 tau=0.1,
                 tauP=0.2,
                 P=1.0,
                 stop_after_percent=100.0,
                 percent_bonds_per_step=0.005,
                 AB_bond_const=100,
                 AB_bond_dist=1.0,
                 shrink_time=1e6,
                 shrinkT=2.0,
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
        self.Anp_interaction = Anp_interaction
        self.Bnp_interaction = Bnp_interaction
        self.Cnp_interaction = Cnp_interaction
        self.npnp_interaction = npnp_interaction
        self.Arc_interaction = Arc_interaction
        self.Brc_interaction = Brc_interaction
        self.Crc_interaction = Crc_interaction
        self.nprc_interaction = nprc_interaction
        self.rcrc_interaction = rcrc_interaction
        self.AA_alpha = AA_alpha
        self.AB_alpha = AB_alpha
        self.AC_alpha = AC_alpha
        self.BC_alpha = BC_alpha
        self.num_spheres = num_spheres
        self.radius_sphere = radius_sphere
        self.radius_a = radius_a
        self.radius_b = radius_b
        self.AB_bond_const=AB_bond_const
        self.AB_bond_dist=AB_bond_dist
        self.stop_after_percent = stop_after_percent
        self.percent_bonds_per_step = percent_bonds_per_step
        self.shrink_time = shrink_time
        self.shrinkT = shrinkT
        self.integrator = cmn.Integrators[integrator]
        self.tau = tau
        self.tauP = tauP
        self.P = P
        self._exclude_bonds_from_nlist = True
        self._all_particle_types = None

    def exclude_bonds_from_nlist(self):
        return self._exclude_bonds_from_nlist

    def get_log_quantities(self):
        log_quantities = super().get_log_quantities()+["pair_lj_energy"]
        return log_quantities

    def get_msd_groups(self):
        self.group_a = hoomd.group.type(name='a-particles', type='A')
        self.group_b = hoomd.group.type(name='b-particles', type='B')
        self.group_c = hoomd.group.type(name='nano-particles', type='np')
        msd_groups = [self.group_a, self.group_b, self.group_c]
        return msd_groups

    def get_non_bonded_neighbourlist(self):
        nl = md.nlist.tree()  # cell()
        nl.reset_exclusions(exclusions=['bond', 'body']);
        return nl

    def get_system_from_file(self, file_path, use_time_step_from_file):
        print("############### get_system_from_file called for",file_path)
        if use_time_step_from_file:
            time_step = None
        else:
            time_step = 0  # start simulation from start.
        if file_path.endswith('.hoomdxml'):
            #system = hoomd.deprecated.init.read_xml(file_path, time_step=time_step)
            system = init_rigid.init_wrapper(xmlfile=file_path,
                                             restart_rigid=True,  # note that this is True and the rigid info written
                                             # earlier is being used to initialize the orientations
                                             rigid_flex_xyz_file="rigid_center_flex.xml",
                                             rigid_json_file="rigid_info.json")
        elif file_path.endswith('.gsd'):
            raise ValueError('Reading the most recent frame from gsd file is not yet implemented!')
            system = hoomd.init.read_gsd(file_path, frame=0, time_step=time_step)
        else:
            raise ValueError('No such file as {} exist on disk!'.format(file_path))
        return system

    def finalize_stage(self, stage):
        super().finalize_stage(stage)
        rigid = hoomd.group.rigid_center()
        nonrigid = hoomd.group.nonrigid()
        deprecated.dump.xml(filename=os.path.join(self.output_dir, "rigid_center_flex.xml"),
                            group=hoomd.group.union("both", rigid, nonrigid),
                            position=True,
                            image=True,
                            mass=True,
                            diameter=True,
                            type=True,
                            body=True,
                            orientation=True,
                            inertia=True)

        
    def set_initial_structure(self):
        print('========INITIAIZING STRUCTURE==========')
        num_beads_in_sphere= 1 + (4*((2*self.radius_sphere)**2))
        Sphere.mass = num_beads_in_sphere
        Sphere.name = 'Sphere'
        mbuild_radius = self.radius_sphere / 10
        
        
        volume_bead = (4/3) * pi * (self.radius_a ** 3)
        volume_particles = (self.num_a * volume_bead) + (self.num_b * volume_bead)
        volume_spheres = self.num_spheres * ((4/3) * pi * ((self.radius_sphere + 0.5) ** 3))
        vol_total = volume_particles + volume_spheres
        vol_fraction_spheres = volume_spheres / vol_total
        
        if self.num_spheres > 0:
        #Volume-fraction based calculation of appropriate box volume
            desired_box_volume = vol_total / vol_fraction_spheres
        else:
            desired_box_volume = ((A.mass*self.num_a) + (B.mass*self.num_b)) / self.density
            
        desired_box_dim = ((desired_box_volume) ** (1./3.)) 
        print('desired box dim =', desired_box_dim)
        reduced_density = self.density / 10
        ex_box_vol = (vol_total) / reduced_density
        expanded_box_dim = (ex_box_vol ** (1. / 3.))                   
        print('expanded box dim =', expanded_box_dim)
        
        #mbuild factor of 10
        half_L = expanded_box_dim/2 
        
        print('half_L', half_L)
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

            mix_box = mb.packing.fill_box([A(), B(), Sphere(radius=mbuild_radius), AB()],
                                          [self.num_a, self.num_b, self.num_spheres, 1],
                                          [self.radius_a, self.radius_b, self.radius_sphere, 0.5],
                                          box=box,
                                          fix_orientation=[False, False, True, False])  # ,overlap=0.5)
            if self.num_spheres > 0: 
                mix_box.label_rigid_bodies(discrete_bodies = 'Sphere')

            if self.init_file_name.endswith('.hoomdxml'):
                mix_box.save(self.init_file_name, overwrite=True)
            elif self.init_file_name.endswith('.gsd'):
                mix_box.save(self.init_file_name, write_ff=False, overwrite=True)

            if self.init_file_name.endswith('.hoomdxml'):
                system = cme_utils.manip.convert_rigid.init_wrapper(xmlfile=self.init_file_name)
            elif self.init_file_name.endswith('.gsd'):
                system = hoomd.init.read_gsd(self.init_file_name)

            print('Initial box dimension: {}'.format(system.box.dimensions))

        self._all_particle_types = system.particles.types

        self.nl = self.get_non_bonded_neighbourlist()
        if self.nl is None:
            raise Exception('Neighbourlist is not set')
      
        #start of mixing stage
        
        self.setup_force_fields(stage=cmn.Stages.MIXING)
        size_variant = variant.linear_interp([(0, system.box.Lx), (self.shrink_time, desired_box_dim)])
        md.integrate.mode_standard(dt=self.mix_dt)
        rigid = hoomd.group.rigid_center()
        nonrigid = hoomd.group.nonrigid()
        both_group = hoomd.group.union("both", rigid, nonrigid)
             
        md.integrate.langevin(group=both_group,
                              kT=self.shrinkT,
                              seed=1223445)
        #integrator.set_gamma('A', gamma=self.gamma)
        #integrator.set_gamma('B', gamma=self.gamma)
        #integrator.set_gamma('C', gamma=self.gamma)

        resize = hoomd.update.box_resize(L=size_variant)
        hoomd.run(self.shrink_time)
        snapshot = system.take_snapshot()
        print('Initial box dimension after shrink: {}'.format(snapshot.box))
        return system

    def setup_force_fields(self, stage):
        if self.DEBUG:
            print('=============force fields parameters==============')
            print('self.AB_bond_const', self.AB_bond_const)
            print('self.AB_bond_dist', self.AB_bond_dist)
            print('self.AA_interaction', self.AA_interaction)
            print('self.AB_interaction', self.AB_interaction)
            print('self.gamma', self.gamma)
        if self.num_b > 0 and self.num_a > 0:
            harmonic = md.bond.harmonic()
            harmonic.bond_coeff.set('C-C', k=self.CC_bond_const, r0=self.CC_bond_dist)
            harmonic.bond_coeff.set('A-B', k=self.AB_bond_const, r0=self.AB_bond_dist)
        lj = md.pair.lj(r_cut=2.5, nlist=self.nl)
        
        lj.pair_coeff.set('A', 'A', epsilon=self.AA_interaction, sigma=1.0, alpha=self.AA_alpha)
        lj.pair_coeff.set('B', 'B', epsilon=self.AA_interaction, sigma=1.0, alpha=self.AA_alpha)
        lj.pair_coeff.set('np', 'np', epsilon=self.npnp_interaction, sigma=1.0, alpha=self.AA_alpha)
        
        lj.pair_coeff.set('A', 'rc', epsilon=0, sigma=0, alpha=self.AA_alpha, r_cut=0)
        lj.pair_coeff.set('B', 'rc', epsilon=0, sigma=0, alpha=self.AA_alpha, r_cut=0)      
        lj.pair_coeff.set('rc', 'rc', epsilon=0, sigma=0, alpha=self.AA_alpha, r_cut=0)
        lj.pair_coeff.set('np', 'rc', epsilon=0, sigma=0, alpha=self.AA_alpha, r_cut=0)
        
        lj.pair_coeff.set('A', 'B', epsilon=self.AB_interaction, sigma=1.0, alpha=self.AB_alpha)
        lj.pair_coeff.set('A', 'np', epsilon=self.Anp_interaction, sigma=1.0, alpha=self.AC_alpha)
        lj.pair_coeff.set('B', 'np', epsilon=self.Bnp_interaction, sigma=1.0, alpha=self.BC_alpha)
        

        # the list comprehension below is setting the interaction parameter for all the particles in the system and
        # '_R' (the rigid body centers) to be zero. As a caveat, don't create a particle type '_R'!
        lj.pair_coeff.set(self._all_particle_types,
                          [i for (i, v) in zip(self._all_particle_types,
                                               [_.startswith("_R") for _ in self._all_particle_types]) if v],
                          epsilon=0.0,
                          sigma=0.0,
                          r_cut=0)

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
        nonrigid = hoomd.group.nonrigid()
        both_group = hoomd.group.union("both", rigid, nonrigid)
        if self.integrator == cmn.Integrators.LANGEVIN:
            integrator = md.integrate.langevin(group=both_group,
                                               kT=temperature,
                                               seed=1223445,
                                               noiseless_t=False,
                                               noiseless_r=False)
            integrator.set_gamma('A', gamma=self.gamma)
            integrator.set_gamma('B', gamma=self.gamma)
        elif self.integrator == cmn.Integrators.NPT:
            raise NotImplementedError('The NPT integrator with ABCNPTypeEpoxyLJHarmonicSimulation is not tested.')
            integrator = md.integrate.npt(group=both_group,
                                          tau=self.tau,
                                          tauP=self.tauP,
                                          P=self.P,
                                          kT=temperature)
