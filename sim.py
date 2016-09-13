import hoomd
import hoomd.deprecated

import numpy as np

hoomd.context.initialize('')


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

        for cord in np.linspace(0.0, END, num = self.N):
            type_dic = {
            "A": [cord, 1.0, 0.0],
            "B": [0.0, cord, 1.0],
            "C": [1.0, 0.0, cord]
            }
            self.cords.append(type_dic[self.btype])


def gen_lattice(basis_list):
    R""" Wraper function for hoomd.lattice.unitcell

    Args:
        num_a (int): Number of A type per unit cell
        num_b (int): Number of B type per unit cell
        num_c (int): Number of C type per unit cell
        N (int): Number of unit cells

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
    print(all_cords)
    uc = hoomd.lattice.unitcell(N = N,
                                a1 = [1,0,0],
                                a2 = [0,1,0],
                                a3 = [0,0,1],
                                dimensions = 3,
                                position = all_cords,
                                type_name = all_types,
                                mass = all_masses,
                                charge = all_charges,
                                diameter = all_diameters
                                )
    return uc





n_cells = 1

a = Basis(N = 2)
b = Basis(btype = "B", N = 1, diameter = 2.3)
c = Basis(btype = "C", N = 1, mass = 1.5)


uc = gen_lattice([a,b,c])

system = hoomd.init.create_lattice(unitcell=uc, n=n_cells);
hoomd.deprecated.dump.xml(group = hoomd.group.all(), filename = "out.hoomdxml", all=True)
