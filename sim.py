import hoomd
import sys
from hoomd import deprecated
from hoomd import dump
from hoomd import md

import numpy as np


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




def init_bonds(indexA, indexB):
    snapshot = system.take_snapshot(bonds=True)
    snapshot.bonds.resize(1)
    snapshot.bonds.group[0] = [indexA, indexB]
    snapshot.bonds.types = ['A-B']
    snapshot.bonds.typeid[0] = 0;
    system.restore_snapshot(snapshot)
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('A-B', k=330.0, r0=1.99)



def make_bond(indexA, indexB):
    snapshot = system.take_snapshot(bonds=True)
    n_bonds = snapshot.bonds.N
    snapshot.bonds.resize(n_bonds + 1)
    snapshot.bonds.group[n_bonds] = [indexA, indexB]
    snapshot.bonds.typeid[0] = 0
    system.restore_snapshot(snapshot)



def my_callback(timestep):
    n_bonds = system.bonds.bdata.getN()
    if n_bonds == 0:
        init_bonds(0, 1)
    else:
        make_bond(n_bonds, n_bonds+1)
        print("call back")



if __name__ == "__main__":


    hoomd.context.initialize()

    n_cells = 10

    a = Basis(N = 1)
    b = Basis(btype = "B", N = 1)
    #c = Basis(btype = "C", N = 5)


    uc = gen_lattice([a,b])

    system = hoomd.init.create_lattice(unitcell=uc, n=n_cells);
    # Better api to manipulate snapshots
    snapshot = system.take_snapshot(bonds=True)
    # Gives us some bonds to play with
    #print(system.bonds)
    #print(snapshot.bonds.N)
    #snapshot.bonds.resize(1000)
    #print(snapshot.bonds.N)
    #system.restore_snapshot(snapshot)
    #print(system.bonds)
    #snapshot = system.take_snapshot()
    #print(snapshot.bonds.N)
    kT = 10.0
    nl = md.nlist.cell()
    dpd = md.pair.dpd(r_cut=3.0, nlist=nl, kT=kT, seed=0)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['A', 'B', 'C'], A=10.0, gamma = 1.2)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['B', 'C', 'A'], A=10.0, gamma = 1.2)
    dpd.pair_coeff.set(['A', 'B', 'C'], ['C', 'A', 'B'], A=10.0, gamma = 1.2)

    #harmonic = md.bond.harmonic()


    md.integrate.mode_standard(dt=0.02)
    md.integrate.nve(group=hoomd.group.all())
    deprecated.dump.xml(group = hoomd.group.all(), filename = "in.hoomdxml", all=True)
    hoomd.analyze.callback(callback = my_callback, period = 1e2)
    hoomd.run(1e3)
    #[print(p) for p in system.particles]
    deprecated.dump.xml(group = hoomd.group.all(), filename = "out.hoomdxml", all=True)
