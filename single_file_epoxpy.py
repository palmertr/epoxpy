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
random.seed(1020)

MAX_A_BONDS = 4
MAX_B_BONDS = 2
rank_dict = {}
get_rank = rank_dict.get
cut_off_dist = 1


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


def print_pos(step):
    snapshot = system.take_snapshot(bonds=True)

    with open("test.txt", "a") as myfile:
        myfile.write('particle 0 position:{} ts:{}'.format(snapshot.particles.position[0], step))

    system.restore_snapshot(snapshot)
    #print(snapshot.bonds.group)

    #print('particle 0 position:{} ts:{}'.format(snapshot.particles.position[0], step))


def make_bond(index_a, index_b, snapshot):
    n_bonds = snapshot.bonds.N
    snapshot.bonds.resize(n_bonds + 1)
    snapshot.bonds.group[n_bonds] = [index_a, index_b]
    # sets new bond to be A-B type
    snapshot.bonds.typeid[n_bonds] = 1  # we know A-B bond type's id is 1


def get_bond_rank(index, snapshot):
    return np.count_nonzero(snapshot.bonds.group == index)
    #rank = 0
    ##print('bonds.group type: {}'.format(type(snapshot.bonds.group)))
    ##print('num bond groups: {}'.format(snapshot.bonds.N))
    #for bond in snapshot.bonds.group:
    #    if bond[0] == index or bond[1] == index:
    #        rank += 1
    #return rank


def pbc_diff(p1, p2, axes):
    dr = p1 - p2
    for i, p in enumerate(dr):
        dr[i] = abs(dr[i])
        if dr[i] > axes[i]*0.5:
            dr[i] -= axes[i]
    return np.sqrt(dr[0]**2+dr[1]**2+dr[2]**2)


def bond_test(kT, delta_e, bond_rank):
    # No idea if this is thread safe, watch out for MPI gotchas
    # Should be able to to tune rate with delta_e and kT
    # delta_e = 1
    mb_stats = np.exp(-delta_e / kT)
    # Divides by bond rank to make it less probable, add one to prevent rank 0
    # issues
    weight = 1
    if bond_rank >= 1:
        weight = 500
    if mb_stats / float(weight) > random.random():
        return True
    else:
        return False

def make_bond(index_a, index_b, snapshot):
    n_bonds = snapshot.bonds.N
    snapshot.bonds.resize(n_bonds + 1)
    snapshot.bonds.group[n_bonds] = [index_a, index_b]
    # sets new bond to be A-B type
    snapshot.bonds.typeid[n_bonds] = 1  # we know A-B bond type's id is 1

def get_bond_rank(index, snapshot):
    return np.count_nonzero(snapshot.bonds.group == index)
    #rank = 0
    ##print('bonds.group type: {}'.format(type(snapshot.bonds.group)))
    ##print('num bond groups: {}'.format(snapshot.bonds.N))
    #for bond in snapshot.bonds.group:
    #    if bond[0] == index or bond[1] == index:
    #        rank += 1
    #return rank

def pbc_diff(p1, p2, axes):
    dr = p1 - p2
    for i, p in enumerate(dr):
        dr[i] = abs(dr[i])
        if dr[i] > axes[i]*0.5:
            dr[i] -= axes[i]
    return np.sqrt(dr[0]**2+dr[1]**2+dr[2]**2)

def bond_test(kT, delta_e, bond_rank):
    # No idea if this is thread safe, watch out for MPI gotchas
    # Should be able to to tune rate with delta_e and kT
    # delta_e = 1
    mb_stats = np.exp(-delta_e / kT)
    # Divides by bond rank to make it less probable, add one to prevent rank 0
    # issues
    weight = 1
    if bond_rank >= 1:
        weight = 500
    if mb_stats / float(weight) > random.random():
        return True
    else:
        return False

def find_neighbours_and_bond(axis, bond_from_idx, bond_from_type, bond_to_group, bond_to_max_rank,
                             bond_to_type, snapshot, xyz0):
    made_bonds = False
    for p in bond_to_group:
        xyz1 = p.position
        r = pbc_diff(np.array(xyz0), np.array(xyz1), axis)
        if r < cut_off_dist:
            bond_to_idx = p.tag
            bond_to_rank = get_rank(bond_to_idx, -1)
            if bond_to_rank < 0:  # -1 is the default value returned when id is not in dict.
                bond_to_rank = get_bond_rank(bond_to_idx, snapshot)
                rank_dict[bond_to_idx] = 0

            if bond_to_rank < bond_to_max_rank:
                delta_e = 0.1
                if bond_test(log.query('temperature'), delta_e, bond_to_rank):
                    # bond_test
                    make_bond(bond_from_idx, bond_to_idx, snapshot)
                    made_bonds = True
                    rank_dict[bond_from_idx] += 1
                    rank_dict[bond_to_idx] += 1
                    break
                    #print("Found one, bonding {} ({}) to {} ({})".format(bond_from_type, bond_from_idx,
                    #                                                     bond_to_type, bond_to_idx))
                    # print("Rank of A {} type of A {}".format(rank, typeA))
                    # print("Rank of B {} type of B {}".format(bond_rank, typeB))

    return made_bonds

def find_pair(timestep):
    # Until I can hack a way to get access to the neighborlist
    snapshot = system.take_snapshot(bonds=True)
    n_p = snapshot.particles.N
    bond_from_type = "C"
    bond_to_group = None
    bond_to_max_rank = None
    bond_from_max_rank = None
    bond_to_type = None
    while bond_from_type == "C":
        # Keep in mind we if N = 5, index 0..4
        bond_from_idx = random.randint(0, n_p - 1)
        # TODO: This can waste time by selecting a 'C' type
        bond_from_typeid = snapshot.particles.typeid[bond_from_idx]
        bond_from_type = snapshot.particles.types[bond_from_typeid]

    if bond_from_type == "A":
        bond_to_group = msd_groups[2]#.group_b
        bond_to_type = 'B'
        bond_to_max_rank = MAX_B_BONDS
        bond_from_max_rank = MAX_A_BONDS
        # typeB = "B"
        # MAX_RANK_B = self.MAX_B_BONDS
    elif bond_from_type == 'B':
        bond_to_group = msd_groups[0]#.group_a
        bond_to_type = 'A'
        bond_to_max_rank = MAX_A_BONDS
        bond_from_max_rank = MAX_B_BONDS
        # typeB = "A"
        # MAX_RANK_B = self.MAX_A_BONDS
    # Check to see if it can make more bonds

    bond_from_rank = get_rank(bond_from_idx, -1)
    if bond_from_rank < 0:  # -1 is the default value returned when id is not in dict.
        bond_from_rank = get_bond_rank(bond_from_idx, snapshot)
        rank_dict[bond_from_idx] = 0

    if bond_from_rank < bond_from_max_rank:
        xyz0 = snapshot.particles.position[bond_from_idx]
        axis = [snapshot.box.Lx, snapshot.box.Ly, snapshot.box.Lz]

    #           start_time = time.time()
        made_bonds = find_neighbours_and_bond(axis, bond_from_idx, bond_from_type, bond_to_group,
                                                   bond_to_max_rank, bond_to_type, snapshot, xyz0)
            #           stop_time = time.time()
        # print("time to calc neighbours for 1 frame = {}".format(stop_time - start_time))
        if made_bonds is True:
            system.restore_snapshot(snapshot)

            #print(snapshot.bonds.group)
            #exit()
    #print('particle 0 position:{} ts:{}'.format(snapshot.particles.position[0], timestep))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if '1' == sys.argv[1]:
            ext_init = True
        elif '2' == sys.argv[1]:
            ext_init = False
            save_bench = False
        else:
            ext_init = False
            save_bench = True
    else:
        ext_init = False
    expected_gsd = 'expected.gsd'
    print("ext_init", ext_init)

    if ext_init:
        ext_init_struct_path = 'initial_bench.hoomdxml' #None
    else:
        ext_init_struct_path = None
    n_mul = 1
    box = [3, 3, 3]
    init_file_name = 'initial.hoomdxml'
    output_dir = '.'
    dcd_write = 1e2
    dt = 1e-2
    mix_time = 3e4
    mix_kT = 2.0
    cure_kt = 2.0
    time_scale = 100
    temp_prof = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kT, initial_time=mix_time)
    temp_prof.add_state_point(500 * time_scale, cure_kt)
    bond = True
    log_write = 1e5
    bond_period = 1e1
    legacy_bonding = True
    log_curing = False
    final_time = temp_prof.get_total_sim_time()
    md__total_time = final_time - mix_time
    md_time = md__total_time
    shrink = False
    num_a = 10 * n_mul
    num_b = 20 * n_mul
    num_c10 = 2 * n_mul
    density = 1.0
    shrink_time = 1.0

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
        log = hoomd.analyze.log(filename=None, quantities=["temperature"], period=bond_period)
        bond_callback = hoomd.analyze.callback(callback=find_pair, period=bond_period)
        #if legacy_bonding is True:
        #    bonding_callback = LegacyBonding(system=system, groups=msd_groups, log=log)
        #else:
        #    bonding_callback = FreudBonding(system=system, groups=msd_groups, log=log)
        #bond_callback = hoomd.analyze.callback(callback=bonding_callback, period=bond_period)

    #    if log_curing is True:
    #        curing_callback = hoomd.analyze.callback(callback=self.calculate_curing_percentage,
    #                                                 period=self.curing_log_period, phase=-1)

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

    if bond is True:
        bond_callback.disable()
        log.disable()
     #   if log_curing is True:
     #       curing_callback.disable()
    print("ext_init", ext_init)
    if ext_init:
        same_same = filecmp.cmp('final.hoomdxml', 'final_bench.hoomdxml')
        assert same_same == True

    if save_bench:
        shutil.copy('/Users/stephenthomas/projects/epoxy_sim/final.hoomdxml',
                    '/Users/stephenthomas/projects/epoxy_sim/final_bench.hoomdxml')
        shutil.copy('/Users/stephenthomas/projects/epoxy_sim/initial.hoomdxml',
                    '/Users/stephenthomas/projects/epoxy_sim/initial_bench.hoomdxml')

    #f = gsd.fl.GSDFile('data.gsd', 'rb')
    #t = gsd.hoomd.HOOMDTrajectory(f)
    #snapshot = t[-1]
    #f = gsd.fl.GSDFile(expected_gsd, 'rb')
    #t = gsd.hoomd.HOOMDTrajectory(f)
    #expected_snapshot = t[-1]
    #assert snapshot.particles.N == expected_snapshot.particles.N
    #expected_pos = expected_snapshot.particles.position
    #current_pos = snapshot.particles.position
    #print('Positions. benchmark:{}, current simulation:{}'.format(expected_pos[0],
    #                                                                                        current_pos[0]))
    #assert np.allclose(expected_pos, current_pos)
    #expected_bonds = expected_snapshot.bonds.N
    #current_bonds = snapshot.bonds.N
    #print('test_epoxy_sim_legacy_bonding_count. benchmark:{}, current simulation:{}'.format(expected_bonds,
    #                                                                                    current_bonds))
    #assert current_bonds == expected_bonds
