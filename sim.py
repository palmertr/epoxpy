import hoomd
import hoomd.deprecated

import numpy as np

hoomd.context.initialize('')


def gen_pos(n_a, n_b, n_c):

    END = 1.0
    a_pos = np.linspace(0.0, END, num = n_a)
    b_pos = np.linspace(0.0, END, num = n_b)
    c_pos = np.linspace(0.0, END, num = n_c)

    some_pos = [a_pos, b_pos, c_pos]

    cords = []
    for axis in (0, 1, 2):
        cord = []
        for pos in some_pos[axis]:
            axis_dic = {
                    0 : [pos, 1.0, 0.0],
                    1 : [0.0, pos, 1.0],
                    2 : [1.0, 0.0, pos]
                    }
            cord.append(axis_dic[axis])
        cords.append([c for c in cord])

    # some flatten magic
    # http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python/952952#952952
    return [item for sublist in cords for item in sublist]

def gen_type_name(n_a, n_b, n_c):
    types = []
    for n in range(n_a):
        types.append("A")
    for n in range(n_b):
        types.append("B")
    for n in range(n_c):
        types.append("C")
    return types


def gen_lattice(num_a, mass_a, dia_a, num_b, mass_b, dia_b, num_c, mass_c, dia_c, N_cells):
    R""" Wraper function for hoomd.lattice.unitcell

    Args:
        n_a (int): Number of A type per unit cell
        n_b (int): Number of B type per unit cell
        n_c (int): Number of C type per unit cell
        N (int): Number of unit cells

    Returns:
        hoomd.lattice.unitcell

    """
    my_cord = gen_pos(n_a, n_b, n_c)
    my_types = gen_type_name(n_a, n_b, n_c)

    print(my_cord)
    print(my_types)
    print([1.0 for i in range(N)])
    uc = hoomd.lattice.unitcell(N = N,
                                a1 = [1,0,0],
                                a2 = [0,1,0],
                                a3 = [0,0,1],
                                dimensions = 3,
                                position = my_cord,
                                type_name = my_types,
                                mass = [mass_a]*num_a,[mass_b]*num_b,[mass_c]*num_c,
                                charge = [0.0]*(num_a+num_b+num_c),
                                diameter = [dia_a]*num_a,[dia_b]*num_b,[dia_c]*num_c
                                )
    return uc




N_cells = 10

num_a = 2
mass_a = 1.0
dia_a = 1.0

num_b = 2
mass_b = 1.0
dia_b = 1.0

num_c = 2
mass_c = 1.0
dia_c = 1.0

gen_lattice(num_a, mass_a, dia_a, num_b, mass_b, dia_b, num_c, mass_c, dia_c, N_cells)

zeros = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
uc = hoomd.lattice.unitcell(N = 3,
                            a1 = [1,0,0],
                            a2 = [0,1,0],
                            a3 = [0,0,1],
                            dimensions = 3,
                            position = [[0,0,0], [0, 0, 1], [1,0,0]],
                            type_name = ['A', 'B', 'C'],
                            mass = [1.0, 1.0, 1.0],
                            charge = [0.0, 0.0, 0.0],
                            diameter = [1.0, 1.0, 1.0])
                            #moment_inertia = zeros,
                            #orientation = zeros);



system = hoomd.init.create_lattice(unitcell=uc, n=10);
hoomd.deprecated.dump.xml(group = hoomd.group.all(), filename = "out.hoomdxml", all=True)

