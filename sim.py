from __future__ import division
import os
import sys
import random
import shutil
import subprocess as sp
import numpy as np
import hoomd
import init as my_init
from hoomd import deprecated
from hoomd import dump
from hoomd import md


def make_bond(indexA, indexB, snapshot):
    n_bonds = snapshot.bonds.N
    snapshot.bonds.resize(n_bonds + 1)
    snapshot.bonds.group[n_bonds] = [indexA, indexB]
    #sets new bond to be A-B type
    snapshot.bonds.typeid[n_bonds] = 1
    system.restore_snapshot(snapshot)


def bond_test(kT, delta_e, bond_rank):
    # No idea if this is thread safe, watch out for MPI gotchas
    # Should be able to to tune rate with delta_e and kT
    #delta_e = 1
    mb_stats = np.exp(-delta_e/kT)
    # Devides by bond rank to make it less probable, add one to prvent rank 0
    # issues
    weight = 1
    if bond_rank >= 1:
        weight = 500
    if mb_stats/float(weight) > random.random():
        return True
    else:
        return False



def get_bond_rank(index, snapshot):
    rank = 0
    for bond in snapshot.bonds.group:
        if bond[0] == index or bond[1] == index:
            rank += 1
    return rank



def pbc_diff(p1, p2, axes):
    dr = p1 - p2
    for i, p in enumerate(dr):
        dr[i] = abs(dr[i])
        if dr[i] > axes[i]*0.5:
            dr[i]-=axes[i]
    #return np.linalg.norm(dr)
    return np.sqrt(dr[0]**2+dr[1]**2+dr[2]**2)



def find_pair(timestep):
    # Until I can hack a way to get access to the neighborlist
    snapshot = system.take_snapshot(bonds=True)
    N_p = snapshot.particles.N
    typeA = "C"
    while typeA == "C":
    # Keep in mind we if N = 5, index 0..4
        indexA = random.randint(0, N_p-1)
    #TODO: This can waste time by selecting a 'C' type
        typeA = system.particles[indexA].type
    if typeA == "A":
        group = groupB
        typeB = "B"
        MAX_RANK_B = MAX_B_BONDS
    else:
        group = groupA
        typeB = "A"
        MAX_RANK_B = MAX_A_BONDS
    # Check to see if it can make more bonds

    rank = get_bond_rank(indexA, snapshot)

    if typeA == "A" and rank == MAX_A_BONDS:
        print("Can't bond it anymore")
    elif typeA == "B" and rank == MAX_B_BONDS:
        print("Can't bond it anymore")

    # If its mixed, it shouldn't be too bad to loop till we find one even if
    # we always start loop from indexA to indexA+1 % N_p
    # also this could be an info loop, yolo!
    else:
        xyz0 = system.particles[indexA].position
        axis = [system.box.Lx, system.box.Ly, system.box.Lz]
        for p in group:
            xyz1 = p.position
            r = pbc_diff(np.array(xyz0), np.array(xyz1), axis)
            if r < CUT:
                indexB = p.tag
                bond_rank = get_bond_rank(indexB, snapshot)
                if bond_rank < MAX_B_BONDS:
                    delta_e = 0.1
                    if bond_test(temp_log.query('temperature'), delta_e, bond_rank):
                        #bond_test
                        make_bond(indexA, indexB, snapshot)
                        print("Found one, bonding {} and {}".format(indexA, indexB))
                        #print("Rank of A {} type of A {}".format(rank, typeA))
                        #print("Rank of B {} type of B {}".format(bond_rank, typeB))

def get_system_mass():
    total_mass = 0
    for p in system.particles:
        total_mass += p.mass
    return total_mass/system.box.get_volume()

def calc_target_V(rho_i, rho_f, V_i):
    return rho_i*V_i/rho_f


if __name__ == "__main__":
    # These will be in an infile somday
    hoomd.context.initialize()
    cwd = os.getcwd()
    #run_name_postfix = sys.argv[1]
    run_name = sys.argv[1]
    run_dir = "/"+run_name
    #print(run_name_postfix)
    print(run_dir)
    #end_eql_kT = float(run_name_postfix)

    #Init System
    A = my_init.Bead()
    B = my_init.Bead(btype="B", mass = 1.0)
    C = my_init.PolyBead(btype="C", mass = 1.0, N = 10)
    # 40 wt C = 2,000
    # 10 wt C = 1,667
    #snap = my_init.init_system({A : 2000, B : 4000, C : 0}, 1)
    N_mult = float(sys.argv[2])
    snap = my_init.init_system({A : int(10000*N_mult), B : int(20000*N_mult), C : int(2000*N_mult)}, 1)
    system = hoomd.init.read_snapshot(snap)

    #Sys Parmas
    log_write = 1e4
    dcd_write = 1e4
    elapsed_time = 0
    # MIX

    time_scale = 7e3
    t_scale = float(sys.argv[2])
    mix_time = 3e6
    mix_kT = 1.0 #hoomd.variant.linear_interp(points = [(0, 5.0), (mix_time, 5.0)])
    elapsed_time += mix_time
    run_type = "A"
    if run_type == "A":
    #NO HOLD
        run_kT = hoomd.variant.linear_interp(points = [(mix_time, 1.0*t_scale), (60*time_scale, 4.5*t_scale), (190*time_scale, 4.5*t_scale), (240*time_scale, 1.0*t_scale)])
    #HOLD
    if run_type == "B":
        run_kT = hoomd.variant.linear_interp(points = [(mix_time, 1.0*t_scale), (50*time_scale, 4.0*t_scale), (120*time_scale, 4.0*t_scale), (140*time_scale, 4.5*t_scale), (270*time_scale, 4.5*t_scale), (320*time_scale, 1.0*t_scale)])

    # BOND
    # Bond cut off
    BOND = False
    bond_end_kT = 1.0
    if BOND == True:
        CUT = 1.0
        MAX_A_BONDS = 4
        MAX_B_BONDS = 2
        bond_period = 1e1
        #bond_time = 1e5 #float(run_name_postfix)
        #bond_end_kT = 1.0 #float(sys.argv[3])
        #bond_kT = float(sys.argv[2]) #hoomd.variant.linear_interp(points = [(elapsed_time, 1.0), (bond_time+elapsed_time, bond_end_kT), (bond_time*2+elapsed_time, bond_end_kT)])
        bond_time = (240*7e3) #bond_time*2
        elapsed_time += bond_time
        print("Number of bonding steps: {}".format(bond_time/bond_period))
    else:
        bond_time = 0

    # EQL
    #end_eql_kT = float(sys.argv[3])
    #eql_kT = float(sys.argv[2])#hoomd.variant.linear_interp(points = [(elapsed_time, bond_end_kT), (eql_time+elapsed_time, end_eql_kT), (eql_time*2+elapsed_time, end_eql_kT)])
    #eql_time = 1e5
    #elapsed_time += eql_time
    ###


    #Mix Step/MD Setup
    groupA = hoomd.group.type(name='a-particles', type='A')
    groupB = hoomd.group.type(name='b-particles', type='B')
    groupC = hoomd.group.type(name='c-particles', type='C')
    hoomd.meta.dump_metadata(filename =  cwd + run_dir + "metadata.json", indent=2)
    deprecated.dump.xml(group = hoomd.group.all(), filename =cwd + run_dir + "start.hoomdxml", all=True)
    hoomd.analyze.log(filename= cwd + run_dir + "out.log", quantities=["pair_dpd_energy","volume","momentum","potential_energy","kinetic_energy","temperature","pressure", "bond_harmonic_energy"], period=log_write, header_prefix='#', overwrite=True)
    dump.dcd(filename=cwd + run_dir +"traj.dcd", period=dcd_write, overwrite=True)
    dump.gsd(filename= cwd + run_dir + "data.gsd", period=dcd_write, group = hoomd.group.all(), overwrite=True, static=['attribute'])

    nl = md.nlist.cell()
    dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=mix_kT, seed=0)
    #dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=1.0, gamma = 1.0)
    dpd.pair_coeff.set('A','A', A=1.0, gamma = 1.0)
    dpd.pair_coeff.set('B','B', A=1.0, gamma = 1.0)
    dpd.pair_coeff.set('C','C', A=1.0, gamma = 1.0)

    dpd.pair_coeff.set('A', 'B', A=10.0, gamma = 1.0)
    dpd.pair_coeff.set('A', 'C', A=10.0, gamma = 1.0)
    dpd.pair_coeff.set('B', 'C', A=10.0, gamma = 1.0)

    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('C-C', k=100.0, r0=1.0)
    harmonic.bond_coeff.set('A-B', k=100.0, r0=1.0)

    md.integrate.mode_standard(dt=1e-2)
    md.integrate.nve(group=hoomd.group.all())
    hoomd.run(mix_time)
    deprecated.dump.xml(group = hoomd.group.all(), filename = cwd + run_dir +"mix.hoomdxml", all=True)
    dpd.set_params(kT = run_kT)
    # Now we bond!
    if BOND is True:
        #dpd.set_params(kT = bond_kT)
        bond_callback = hoomd.analyze.callback(callback = find_pair, period = bond_period)
        temp_log = hoomd.analyze.log(filename=None, quantities=["temperature"], period = bond_period)
        hoomd.run(bond_time)
        bond_callback.disable()
        temp_log.disable()
        deprecated.dump.xml(group = hoomd.group.all(), filename =cwd +run_dir + "bond.hoomdxml", all=True)
    # Now we run to eql
    #dpd.set_params(kT = eql_kT)
    deprecated.analyze.msd(groups=[groupA, groupB, groupC], period=log_write, filename= cwd + run_dir + "msd.log", header_prefix='#')
    #NO HOLD
    if run_type == "A":
        hoomd.run(2e6)
    #HOLD
    if run_type == "B":
        hoomd.run(2e6)
    deprecated.dump.xml(group = hoomd.group.all(), filename = cwd +run_dir +"final.hoomdxml", all=True)
    print("sim fin")
