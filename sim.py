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
    print(snapshot.bonds.types)
    #snapshot.bonds.types = ['A-B'] #Shouldn't do this every time
    #sets new bond to be A-B type
    snapshot.bonds.typeid[n_bonds] = 1
    system.restore_snapshot(snapshot)
    # TODO: This overwrites C-C bond information


def bond_test(kT, delta_e, bond_rank):
    # No idea if this is thread safe, watch out for MPI gotchas
    # Should be able to to tune rate with delta_e and kT
    #delta_e = 1
    mb_stats = np.exp(-delta_e/kT)
    # Devides by bond rank to make it less probable, add one to prvent rank 0
    # issues
    if mb_stats/float(bond_rank+1) > random.random():
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
    return np.linalg.norm(dr)


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
        return False
    elif typeA == "B" and rank == MAX_B_BONDS:
        print("Can't bond it anymore")
        return False

    # If its mixed, it shouldn't be too bad to loop till we find one even if
    # we always start loop from indexA to indexA+1 % N_p
    # also this could be an info loop, yolo!
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
                kT = bond_kT
                if bond_test(kT, delta_e, bond_rank):
                    #bond_test
                    make_bond(indexA, indexB, snapshot)
                    print("Found one, bonding {} and {}".format(indexA, indexB))
                    #print("Rank of A {} type of A {}".format(rank, typeA))
                    #print("Rank of B {} type of B {}".format(bond_rank, typeB))
                    return True

def get_system_mass():
    total_mass = 0
    for p in system.particles:
        total_mass += p.mass
    return total_mass/system.box.get_volume()

def calc_target_V(rho_i, rho_f, V_i):
    return rho_i*V_i/rho_f


if __name__ == "__main__":
    # These will be in an infile somday
    cwd = os.getcwd()
    # Do not change this sys.argv!
    run_name_postfix = sys.argv[1]
    run_name = sys.argv[2]
    run_dir = "/"+run_name
    print(run_name_postfix)
    print(run_dir)


    MAX_A_BONDS = 4
    MAX_B_BONDS = 2

    hoomd.context.initialize()
    BOND = True
    CUT = 1.0
    kT = float(run_name_postfix)
    mix_time = 1e4
    mix_kT = 10.0
    bond_kT = 10.0
    log_write = 1e4
    dcd_write = 1e4
    bond_period = 1e1
    bond_time = 1e3
    final_run_time = 1e3
    run_kT = kT

    A = my_init.Bead()
    B = my_init.Bead(btype="B", mass = 1.0)
    C = my_init.PolyBead(btype="C", mass = 1.0, N = 10)
    # 40 wt C = 2,000
    # 10 wt C = 1,667

    snap = my_init.init_system({A : 10, B : 20, C : 2}, 1)

    system = hoomd.init.read_snapshot(snap)

    # Could make this dynamic by getting types first

    groupA = hoomd.group.type(name='a-particles', type='A')
    groupB = hoomd.group.type(name='b-particles', type='B')

    print(cwd + run_dir + "start.hoomdxml")
    deprecated.dump.xml(group = hoomd.group.all(), filename =cwd + run_dir + "start.hoomdxml", all=True)
    hoomd.analyze.log(filename= cwd + run_dir + "out.log", quantities=["pair_dpd_energy","volume","momentum","potential_energy","kinetic_energy","temperature","pressure", "bond_harmonic_energy"], period=log_write, header_prefix='#', overwrite=True)
    dump.dcd(filename=cwd + run_dir +"traj.dcd", period=dcd_write, overwrite=True)
    # Now we need a mix step

    nl = md.nlist.cell()
    dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=mix_kT, seed=0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=1.0, gamma = 1.0)

    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('C-C', k=100.0, r0=1.0)
    harmonic.bond_coeff.set('A-B', k=100.0, r0=1.0)
    
    md.integrate.mode_standard(dt=1e-2)
    md.integrate.nve(group=hoomd.group.all())
    # particles mix
    hoomd.run(mix_time)
    # Set things to normal coefs

    dpd.pair_coeff.set('A','A', A=1.0, gamma = 1.0)
    dpd.pair_coeff.set('B','B', A=1.0, gamma = 1.0)
    dpd.pair_coeff.set('C','C', A=1.0, gamma = 1.0)

    dpd.pair_coeff.set('A', 'B', A=10.0, gamma = 1.0)
    dpd.pair_coeff.set('A', 'C', A=10.0, gamma = 1.0)
    dpd.pair_coeff.set('B', 'C', A=10.0, gamma = 1.0)



    deprecated.dump.xml(group = hoomd.group.all(), filename = cwd + run_dir +"mix.hoomdxml", all=True)

    # Now we bond!
    if BOND is True:
        dpd.set_params(kT = bond_kT)
        bond_callback = hoomd.analyze.callback(callback = find_pair, period = bond_period)
        hoomd.run(bond_time)
        bond_callback.disable()
        deprecated.dump.xml(group = hoomd.group.all(), filename =cwd +run_dir + "bond.hoomdxml", all=True)
    # Now we run to eql
    
    dpd.set_params(kT = run_kT)
    hoomd.run(final_run_time)
    deprecated.dump.xml(group = hoomd.group.all(), filename = cwd +run_dir +"final.hoomdxml", all=True)
    print("sim fin")
