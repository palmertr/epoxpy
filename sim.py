from __future__ import division
import os
import sys
import random
import shutil
import subprocess as sp
import numpy as np
import hoomd
from hoomd import deprecated
from hoomd import dump
from hoomd import md


class Basis:

    def __init__(self, diameter=1.0, mass=1.0, charge=0.0, btype="A", N=1):

        self.diameter = diameter
        self.mass = mass
        self.charge = charge
        self.btype = btype
        self.N = N
        self.gen_cords()


    def gen_cords(self):

        END = 1.0
        self.cords = []

        for cord in np.linspace(0.1, END, num = self.N):
            type_dic = {
            "A": [cord, 0.5, 0.5],
            "B": [0.5, cord, 0.5],
            "C": [0.5, 0.5, cord]
            }
            self.cords.append(type_dic[self.btype])

    def get_basis_mass(self):
        return self.mass*self.N


def gen_lattice(basis_list, rho):
    R""" Wraper function for hoomd.lattice.unitcell

    Args:
        basis_list (list): List of basis to use in unit cell

    Returns:
        hoomd.lattice.unitcell

    """
    all_types = []
    all_cords = []
    all_masses = []
    all_charges = []
    all_diameters = []
    total_mass = 0
    N = 0

    for basis in basis_list:
        all_types += [basis.btype]*basis.N
        all_cords += basis.cords
        all_masses += [basis.mass]*basis.N
        all_charges += [basis.charge]*basis.N
        all_diameters += [basis.diameter]*basis.N
        N += basis.N
        total_mass += basis.get_basis_mass()

    # volume calc

    V = total_mass/rho
    print("total mass is {}".format(total_mass))

    L = V**(1/3)

    uc = hoomd.lattice.unitcell(N = N,
                                a1 = [L, 0, 0],
                                a2 = [0, L, 0],
                                a3 = [0, 0, L],
                                dimensions = 3,
                                position = all_cords,
                                type_name = all_types,
                                mass = all_masses,
                                charge = all_charges,
                                diameter = all_diameters
                                )
    return uc


#TODO Right now every function and their mom makes snapshots, should be able
# to make one and pass it arround

def init_bonds():
    find_pair()
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('A-B', k=330.0, r0=1.99)



def make_bond(indexA, indexB):
    snapshot = system.take_snapshot(bonds=True)
    n_bonds = snapshot.bonds.N
    snapshot.bonds.resize(n_bonds + 1)
    snapshot.bonds.group[n_bonds] = [indexA, indexB]
    snapshot.bonds.types = ['A-B'] #Shouldn't do this every time
    snapshot.bonds.typeid[0] = 0
    system.restore_snapshot(snapshot)
    # TODO: dynamicly set the bond lenght


def bond_test(kT):
    # No idea if this is thread safe, watch out for MPI gotchas
    #random.seed(a=0)
    # Should be able to to tune rate with delta_e and kT
    delta_e = 1
    mb_stats = np.exp(-delta_e/kT)
    if mb_stats > random.random():
        return True
    else:
        return False

def get_distance(xyz0, xyz1):
    return np.linalg.norm(np.array(xyz0)-np.array(xyz1))


def get_bond_rank(index):
    rank = 0
    snapshot = system.take_snapshot(bonds=True)
    for bond in snapshot.bonds.group:
        if bond[0] == index or bond[1] == index:
            rank += 1
    return rank


def find_pair(timestep):
    # Until I can hack a way to get access to the neighborlist
    snapshot = system.take_snapshot()
    N_p = snapshot.particles.N
    # Keep in mind we if N = 5, index 0..4
    indexA = random.randint(0, N_p-1)
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

    rank = get_bond_rank(indexA)

    if typeA == "A" and rank == MAX_A_BONDS:
        print("Can't bond it anymore")
        return False
    elif typeA == "B" and rank == MAX_B_BONDS:
        print("Can't bond it anymore")
        return False

    # If its mixed, it shouldn't be too bad to loop till we find one even if
    # we always start loop from indexA to indexA+1 % N_p
    # also this could be an info loop, yolo!
    found = False
    #print("Going into the loop")
    xyz0 = system.particles[indexA].position
    # Index magic, makes sure we loop over all particles
    for p in group:
        xyz1 = p.position
        # Ugh, need to account for PBC
        r = get_distance(xyz0, xyz1)
        if r < CUT:
            bond_rank = get_bond_rank(p.tag)
            if bond_rank < MAX_B_BONDS:
                #bond_test
                indexB = p.tag
                make_bond(indexA, indexB)
                #print("Found one, bonding {} and {}".format(indexA, indexB))
                found = True
                #print("Rank of A {} type of A {}".format(rank, typeA))
                #print("Rank of B {} type of B {}".format(get_bond_rank(p.tag), typeB))
                return True





def my_callback(timestep):
    n_bonds = system.bonds.bdata.getN()
    if n_bonds == 0:
        init_bonds()
        print("first bond")
    else:
        print("next bond")
        if find_pair() is True:
            print("bonded")
        else:
            print("could not bond")


def get_system_mass():
    total_mass = 0
    for p in system.particles:
        total_mass += p.mass
    return total_mass/system.box.get_volume()

def calc_target_V(rho_i, rho_f, V_i):
    return rho_i*V_i/rho_f

def init_run_dir(run_dir):
    print("making dir..")
    cwd = os.getcwd()
    if os.path.isdir(cwd+ run_dir) is not True:
        print(cwd + run_dir)
        print(os.path.isfile(cwd+ run_dir))
        os.makedirs(cwd + run_dir)
        print("made")
        print(cwd + run_dir)



if __name__ == "__main__":
    # These will be in an infile somday
    cwd = os.getcwd()
    # Do not change this sys.argv!
    run_name_postfix = sys.argv[1]
    run_name = sys.argv[2]
    #run_name = "dpdc_debug_bonding_p10g0_{}/".format(run_name_postfix)
    run_dir = "/"+run_name
    #init_run_dir(run_dir)
    #shutil.copy(cwd + "/submit.sh", cwd + run_dir + "submit.sh")
    #shutil.copy(cwd + "/sim.py", cwd + run_dir + "sim.py")
    print(run_name_postfix)
    print(run_dir)

    # run vars below

    MAX_A_BONDS = 4
    MAX_B_BONDS = 2

    hoomd.context.initialize()
    BOND =False
    CUT = 2.0
    kT = float(run_name_postfix)
    n_cells = 30 #2*30^3 = 54k
    a = Basis(N = 1)
    b = Basis(btype = "B", N = 2)
    #c = Basis(btype = "C", N = 5)
    rho = 1.0 #float(run_name_postfix)
    uc = gen_lattice([a,b], rho)
    mix_time = 1e4
    mix_kT = 10.0
    bond_kT = kT
    log_write = 1e4
    dcd_write = 1e4
    bond_period = 1e1
    bond_time = 5e3
    final_run_time = 5e5
    run_kT = kT
    # Maybe the infile returns a snapshot?
    system = hoomd.init.create_lattice(unitcell=uc, n=n_cells);

    # Could make this dynamic by getting types first

    groupA = hoomd.group.type(name='a-particles', type='A')
    groupB = hoomd.group.type(name='b-particles', type='B')

    print(cwd + run_dir + "start.hoomdxml")
    deprecated.dump.xml(group = hoomd.group.all(), filename =cwd + run_dir + "start.hoomdxml", all=True)
    hoomd.analyze.log(filename= cwd + run_dir + "out.log", quantities=["pair_dpd_energy","volume","momentum","potential_energy","kinetic_energy","temperature","pressure", "bond_harmonic_energy"], period=log_write, header_prefix='#', overwrite=True)
    dump.dcd(filename=cwd + run_dir +"traj.dcd", period=dcd_write, overwrite=True)
    # Now we need a mix step

    nl = md.nlist.cell()
    dpd = md.pair.dpd(r_cut=2.0, nlist=nl, kT=mix_kT, seed=0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['A', 'B', 'C'], A=1.0, gamma = 1.0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['B', 'C', 'A'], A=10.0, gamma = 1.0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=10.0, gamma = 1.0)

    # Test to see if this will fix bondsbeing calculated
    # Manualy bond 2 to create bonds
    # TODO fix this hack
    make_bond(0,1)
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('A-B', k=400.0, r0=1.0)


    md.integrate.mode_standard(dt=2e-2)
    md.integrate.nve(group=hoomd.group.all())
    # For the mix we will set r_cut to be larger, this lets all of our
    # particles mix
    hoomd.run(mix_time)
    # Set r_cut back to what it should be
    dpd.pair_coeff.set(['A', 'B', 'C'], ['A', 'B', 'C'], A=1.0, gamma = 1.0, r_cut = 1.0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['B', 'C', 'A'], A=10.0, gamma = 1.0, r_cut = 1.0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=10.0, gamma = 1.0, r_cut = 1.0)

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
    # Clean up now
    # Move run files to dir
    # I think should be done in bash, then no race conditions

    print("clean sim.py exit")
