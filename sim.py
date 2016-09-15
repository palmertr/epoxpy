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


def gen_lattice(basis_list):
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
    N = 0

    for basis in basis_list:
        all_types += [basis.btype]*basis.N
        all_cords += basis.cords
        all_masses += [basis.mass]*basis.N
        all_charges += [basis.charge]*basis.N
        all_diameters += [basis.diameter]*basis.N
        N += basis.N

    uc = hoomd.lattice.unitcell(N = N,
                                a1 = [1, 0, 0],
                                a2 = [0, 1, 0],
                                a3 = [0, 0, 1],
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


def find_pair():
    # Until I can hack a way to get access to the neighborlist
    snapshot = system.take_snapshot()
    N_p = snapshot.particles.N
    # Keep in mind we if N = 5, index 0..4
    indexA = random.randint(0, N_p-1)
    typeA = system.particles[indexA].type
    if typeA == "A":
        group = groupB
    else:
        group = groupA
    # If its mixed, it shouldn't be too bad to loop till we find one even if
    # we always start loop from indexA to indexA+1 % N_p
    # also this could be an info loop, yolo!
    found = False
    print("Going into the loop")
    while found is False:
        xyz0 = system.particles[indexA].position
        # Index magic, makes sure we loop over all particles
        for p in group:
            xyz1 = p.position
            r = get_distance(xyz0, xyz1)
            if r < CUT:
                #bond_test
                indexB = p.tag
                make_bond(indexA, indexB)
                print("Found one, bonding {} and {}".format(indexA, indexB))
                found = True
                break





def my_callback(timestep):
    n_bonds = system.bonds.bdata.getN()
    if n_bonds == 0:
        init_bonds()
        print("first bond")
    else:
        print("next bond")
        find_pair()


def calc_rho():
    total_mass = 0
    for p in system.particles:
        total_mass += p.mass
    return total_mass/system.box.get_volume()

def calc_target_V(rho_i, rho_f, V_i):
    return rho_i*V_i/rho_f

def init_run_dir(run_dir):
    cwd = os.getcwd()
    if os.path.isfile(run_dir):
        yield
    else:     
        os.makedirs(run_dir)

if __name__ == "__main__":
    # These will be in an infile somday
    hoomd.context.initialize()
    run_dir = "runs/"
    run_name = "test/"
    run_dir += run_name
    init_run_dir(run_dir)
    #print(run_dir)
    # Move run files to dir
    shutil.copy("submit.sh", run_dir)
    shutil.copy("sim.py", run_dir)
    # run vars below
    CUT = 5.0
    kT = 10.0
    n_cells = 10
    a = Basis(N = 1)
    b = Basis(btype = "B", N = 1)
    #c = Basis(btype = "C", N = 5)
    uc = gen_lattice([a,b])
    mix_time = 1e3
    mix_kT = 10.0
    rho = 1.0
    shrink_time = 1e3
    shrink_kT = 10.0
    bond_time = 1e3
    log_write = 1e2
    dcd_write = 1e2
    bond_period = 1e2
    bond_time = 1e3
    final_run_time = 1e1
    # Maybe the infile returns a snapshot?
    system = hoomd.init.create_lattice(unitcell=uc, n=n_cells);

    # Could make this dynamic by getting types first

    groupA = hoomd.group.type(name='a-particles', type='A')
    groupB = hoomd.group.type(name='b-particles', type='B')


    deprecated.dump.xml(group = hoomd.group.all(), filename = "start.hoomdxml", all=True)
    hoomd.analyze.log(filename= run_dir + "out.log", quantities=["pair_dpd_energy","volume","momentum","potential_energy","kinetic_energy","temperature","pressure", "bond_harmonic_energy"], period=log_write, header_prefix='#', overwrite=True)
    dump.dcd(filename=run_dir +"traj.dcd", period=dcd_write, overwrite=True)
    # Now we need a mix step

    nl = md.nlist.cell()
    dpd = md.pair.dpd(r_cut=3.0, nlist=nl, kT=mix_kT, seed=0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['A', 'B', 'C'], A=10.0, gamma = 1.2)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['B', 'C', 'A'], A=10.0, gamma = 1.2)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=10.0, gamma = 1.2)

    # Test to see if this will fix bondsbeing calculated
    # Manualy bond 2 to create bonds
    # TODO fix this hack
    make_bond(0,1)
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('A-B', k=330.0, r0=1.99)


    md.integrate.mode_standard(dt=0.02)
    md.integrate.nve(group=hoomd.group.all())
    hoomd.run(mix_time)
    deprecated.dump.xml(group = hoomd.group.all(), filename = run_dir +"mix.hoomdxml", all=True)

    # Some step to get the correct density

    rho_i = calc_rho()
    V_i = system.box.get_volume()
    start_L = V_i**(1/3)
    target_v = calc_target_V(rho_i, rho, V_i)
    target_L = target_v**(1/3)
    print("Starting at {} shrinking to {}".format(start_L, target_L))
    hoomd.update.box_resize(L = hoomd.variant.linear_interp([(0, start_L), (shrink_time, target_L)]))
    hoomd.run(shrink_time)
    deprecated.dump.xml(group = hoomd.group.all(), filename = run_dir +"shrink.hoomdxml", all=True)


    # Now we bond!

    #bond_callback = hoomd.analyze.callback(callback = my_callback, period = bond_period)
    hoomd.run(bond_time)
    #bond_callback.disable()
    deprecated.dump.xml(group = hoomd.group.all(), filename =run_dir + "bond.hoomdxml", all=True)
    hoomd.run(final_run_time)
    deprecated.dump.xml(group = hoomd.group.all(), filename = run_dir +"final.hoomdxml", all=True)
    # Clean up now
    shutil.copy("job.o", run_dir)
    shutil.copy("job_a.o", run_dir)
    shutil.copy("job_b.o", run_dir)
