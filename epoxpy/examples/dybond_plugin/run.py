from hoomd import dump
from epoxpy.lib import A, B, C, C10, Epoxy_A_10_B_20_C10_2_Blend
import hoomd
from hoomd import md
from hoomd import deprecated
import mbuild as mb
import os
import epoxpy.temperature_profile_builder as tpb
import gsd.hoomd
import numpy as np
import random
import sys
import filecmp
import shutil
import hoomd.dybond_plugin as db
random.seed(1020)

MAX_A_BONDS = 4
MAX_B_BONDS = 2
rank_dict = {}
get_rank = rank_dict.get
cut_off_dist = 1
n_mul = 1000
shrink = True
density = 1.0
shrink_time = 1.0
box = [3, 3, 3]#if shrink is True, this box gets changed according to the desired density
init_file_name = 'initial.hoomdxml'
output_dir = '.'
dcd_write = 1e2
dt = 1e-2
mix_time = 3e4
mix_kT = 2.0
cure_kt = 2.0
time_scale = 1e4
bond = True
log_write = 1e5
bond_period = 1e1
log_curing = True
legacy_bonding = False


def initialize_context():
    try:
        __IPYTHON__
        run_from_ipython = True
    except NameError:
        run_from_ipython = False
    if run_from_ipython:
        print('Initializing HOOMD in ipython')
        hoomd.context.initialize('--mode=cpu')
    else:
        hoomd.context.initialize()



if __name__ == '__main__':
    ext_init = False
    expected_gsd = 'expected.gsd'
    print("ext_init", ext_init)

    if ext_init:
        ext_init_struct_path = 'initial_bench.hoomdxml' #None
    else:
        ext_init_struct_path = None
    temp_prof = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kT, initial_time=mix_time)
    temp_prof.add_state_point(500 * time_scale, cure_kt)
    
    final_time = temp_prof.get_total_sim_time()
    md__total_time = final_time - mix_time
    md_time = md__total_time
    num_a = 10 * n_mul
    num_b = 20 * n_mul
    num_c10 = 2 * n_mul
    # initialize hoomd context
    initialize_context()
    #creating initial structure. There are two possible ways of doing it.
    # 1. Use mbuild to create the initial structure and use hoomd to shrink the volume to desired density
    # 2. Load initial structure created externally.
    if ext_init_struct_path is None:
        # initialize structure
        if shrink is True:
            desired_box_volume = ((A.mass * num_a) + (B.mass * num_b) + (C10.mass * num_c10)) / density
            desired_box_dim = (desired_box_volume ** (1. / 3.))
            box = [desired_box_dim, desired_box_dim, desired_box_dim]
            print('Packing {} A particles ..'.format(num_a))
            mix_box = mb.packing.fill_box(A(), num_a, box=box)
            mix_box = mb.packing.solvate(mix_box, B(), num_b, box=box)
            print('Packing {} B particles ..'.format(num_b))
            mix_box = mb.packing.solvate(mix_box, C10(), num_c10, box=box)
            print('Packing {} C10 particles ..'.format(num_c10))
        else:
            blend = Epoxy_A_10_B_20_C10_2_Blend()
            mix_box = mb.packing.fill_box(blend, n_mul, box=box, overlap=0.050)
        mix_box.save(init_file_name)  # save file to disk so that structure can be loaded into hoomd
        system = hoomd.deprecated.init.read_xml(init_file_name)  # initialize hoomd system
        snapshot = system.take_snapshot(bonds=True)  # take snapshot to change particle type ids and bond ids
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
        system.restore_snapshot(snapshot)

        if shrink is True:
            hoomd.update.box_resize(period=1, L=desired_box_dim)
            hoomd.run(shrink_time)
            snapshot = system.take_snapshot()
            print('Initial box dimension: {}'.format(snapshot.box))

        deprecated.dump.xml(group=hoomd.group.all(), filename=init_file_name, all=True)

        #del system  # needed for re initializing hoomd after randomize
        #initialize_context()
        #system = hoomd.deprecated.init.read_xml(init_file_name)
        #snapshot = system.take_snapshot(bonds=True)
        #snapshot.bonds.types = ['C-C', 'A-B']
        #system.restore_snapshot(snapshot)
    else:
        system = hoomd.deprecated.init.read_xml(ext_init_struct_path, time_step=0)
        snapshot = system.take_snapshot(bonds=True)
        snapshot.bonds.types = ['C-C', 'A-B']
        system.restore_snapshot(snapshot)
        deprecated.dump.xml(group=hoomd.group.all(), filename=init_file_name, all=True)

    del system  # needed for re initializing hoomd after randomize
    initialize_context()
    system = hoomd.deprecated.init.read_xml(init_file_name, time_step=0)
    snapshot = system.take_snapshot(bonds=True)
    snapshot.bonds.types = ['C-C', 'A-B']
    system.restore_snapshot(snapshot)

    #mix
    group_a = hoomd.group.type(name='a-particles', type='A')
    group_b = hoomd.group.type(name='b-particles', type='B')
    group_c = hoomd.group.type(name='c-particles', type='C')
    msd_groups = [group_a, group_b, group_c]
    nl = md.nlist.cell()
    dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=mix_kT, seed=1020)
    dpd.pair_coeff.set('A', 'A', A=1.0, gamma=1.0)
    dpd.pair_coeff.set('B', 'B', A=1.0, gamma=1.0)
    dpd.pair_coeff.set('C', 'C', A=1.0, gamma=1.0)
    dpd.pair_coeff.set('A', 'B', A=10.0, gamma=1.0)
    dpd.pair_coeff.set('A', 'C', A=10.0, gamma=1.0)
    dpd.pair_coeff.set('B', 'C', A=10.0, gamma=1.0)
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('C-C', k=100.0, r0=1.0)
    harmonic.bond_coeff.set('A-B', k=100.0, r0=1.0)

    #configure outputs
    hoomd.meta.dump_metadata(filename=os.path.join(output_dir,
                                                           'metadata.json'), indent=2)
    deprecated.dump.xml(group=hoomd.group.all(),
                        filename=os.path.join(output_dir,
                                              'start.hoomdxml'), all=True)
    hoomd.analyze.log(filename=os.path.join(output_dir, 'out.log'),
                      quantities=["pair_dpd_energy", "volume", "momentum", "potential_energy", "kinetic_energy",
                                  "temperature", "pressure", "bond_harmonic_energy"], period=log_write,
                      header_prefix='#', overwrite=True)
    dump.dcd(filename=os.path.join(output_dir, 'traj.dcd'), period=dcd_write, overwrite=True)
    dump.gsd(filename=os.path.join(output_dir, 'data.gsd'), period=dcd_write, group=hoomd.group.all(),
             overwrite=True, static=['attribute'])

    #run mixing
    md.integrate.mode_standard(dt=dt)
    md.integrate.nve(group=hoomd.group.all())
    hoomd.run(mix_time)
    deprecated.dump.xml(group=hoomd.group.all(), filename='mixed.hoomdxml', all=True)

    del system  # needed for re initializing hoomd after randomize
    initialize_context()
    system = hoomd.deprecated.init.read_xml('mixed.hoomdxml', time_step=0)
    snapshot = system.take_snapshot(bonds=True)
    snapshot.bonds.types = ['C-C', 'A-B']
    system.restore_snapshot(snapshot)

    ##### RUN #####
    group_a = hoomd.group.type(name='a-particles', type='A')
    group_b = hoomd.group.type(name='b-particles', type='B')
    group_c = hoomd.group.type(name='c-particles', type='C')
    msd_groups = [group_a, group_b, group_c]
    nl = md.nlist.cell()
    profile = temp_prof.get_profile()
    print('temperature profile {}'.format(profile.points))
    #dpd.set_params(kT=profile)
    dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=profile, seed=12345)

    dpd.pair_coeff.set('A', 'A', A=1.0, gamma=1.0)
    dpd.pair_coeff.set('B', 'B', A=1.0, gamma=1.0)
    dpd.pair_coeff.set('C', 'C', A=1.0, gamma=1.0)
    dpd.pair_coeff.set('A', 'B', A=10.0, gamma=1.0)
    dpd.pair_coeff.set('A', 'C', A=10.0, gamma=1.0)
    dpd.pair_coeff.set('B', 'C', A=10.0, gamma=1.0)
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('C-C', k=100.0,r0=1.0)
    harmonic.bond_coeff.set('A-B', k=100.0, r0=1.0)

    deprecated.analyze.msd(groups=msd_groups, period=log_write, overwrite=True,
                           filename=os.path.join(output_dir, 'msd.log'), header_prefix='#')

    if bond is True:
        #setup the dynamic bond updater and reuse the neighbour list
        updater = db.update.dybond(nl, group=hoomd.group.all(), period=bond_period)
        updater.set_params(bond_type='A-B',A='A',A_fun_groups=MAX_A_BONDS,B='B',B_fun_groups=MAX_B_BONDS,rcut=1.0,Ea=1.0,alpha=2.0)

    #configure outputs
    hoomd.meta.dump_metadata(filename=os.path.join(output_dir,
                                                           'metadata.json'), indent=2)
    deprecated.dump.xml(group=hoomd.group.all(),
                        filename=os.path.join(output_dir,
                                              'start.hoomdxml'), all=True)
    hoomd.analyze.log(filename=os.path.join(output_dir, 'out.log'),
                      quantities=["pair_dpd_energy", "volume", "momentum", "potential_energy", "kinetic_energy",
                                  "temperature", "pressure", "bond_harmonic_energy"], period=log_write,
                      header_prefix='#', overwrite=True)
    dump.dcd(filename=os.path.join(output_dir, 'traj.dcd'), period=dcd_write, overwrite=True)
    dump.gsd(filename=os.path.join(output_dir, 'data.gsd'), period=dcd_write, group=hoomd.group.all(),
             overwrite=True, static=['attribute'])

    hoomd.util.quiet_status()

    md.integrate.mode_standard(dt=dt)
    md.integrate.nve(group=hoomd.group.all())
    deprecated.dump.xml(group=hoomd.group.all(), filename='justbeforerun.hoomxml', all=True)

    print('md run time: ', md_time)
    hoomd.run(md_time)

    deprecated.dump.xml(group=hoomd.group.all(),
                        filename=os.path.join(output_dir,
                                              'final.hoomdxml'), all=True)


