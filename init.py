import numpy as np

from hoomd import * # I really don't like global imports like this

class Bead:

    def __init__(self, diameter=1.0, mass=1.0, charge=0.0, btype="A", N=1):
        self.N = N
        self.diameter = diameter
        self.mass = mass
        self.charge = charge
        self.btype = btype

class PolyBead(Bead):

    def __init__(self, diameter=1.0, mass=1.0, charge=0.0, btype="C", N=5):
        Bead.__init__(self, diameter, mass, charge, btype, N)
        self.bond_type = "{}-{}".format(self.btype, self.btype)



def init_system(building_block_dic, rho):
    """ building_block_dic should be a dctionay of beads/polybeads : # """

    context.initialize("")
    N = 0
    system_mass = 0

    particle_types = []
    particle_bond_types = []

    typeid_list = []
    mass_list = []
    diameter_list = []
    i = 0
    for bead, N_b in building_block_dic.iteritems():
        N += N_b * bead.N
        particle_types.append(bead.btype)
        system_mass += bead.mass * N_b * bead.N
        mass_list += [bead.mass] * N_b * bead.N
        diameter_list += [bead.diameter] * N_b* bead.N
        typeid_list += [i] * N_b * bead.N
        i += 1

        if bead.N > 1:
            particle_bond_types.append(bead.bond_type)


    L = (system_mass/rho)**(1.0/3.0)


    snap = data.make_snapshot(N=N, box=data.boxdim(L=L), particle_types=particle_types)
    snap.particles.typeid[:] = typeid_list
    snap.particles.mass[:] = mass_list


    return snap

A = Bead()
B = PolyBead(btype="B", mass = 2.0, N = 2)

snap = init_system({A : 4, B : 4}, 1)
print(snap.particles.position)
print(snap.particles.mass)
print(snap.particles.types)
print(snap.particles.typeid[0])
print(snap.particles.typeid[-1])
print(snap.box)
