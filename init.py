import numpy as np
import random
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
        self.bond_list = np.array([_ for _ in range(N)])



def init_system(building_block_dic, rho):
    """ building_block_dic should be a dctionay of beads/polybeads : # """

    context.initialize("")
    N = 0
    system_mass = 0

    particle_types = []
    bond_types = []

    typeid_list = []
    mass_list = []
    diameter_list = []
    point_list = []
    bond_list = []
    i = 0
    offset = 0
    for bead, N_b in building_block_dic.iteritems():
        N += N_b * bead.N
        particle_types.append(bead.btype)
        system_mass += bead.mass * N_b * bead.N
        mass_list += [bead.mass] * N_b * bead.N
        diameter_list += [bead.diameter] * N_b* bead.N
        typeid_list += [i] * N_b * bead.N


        if bead.N > 1:
            bond_types.append(bead.bond_type)
            for k in range(N_b):
                bnd_list_str = bead.bond_list +  k*bead.N + offset
                for a,b in zip(bnd_list_str, bnd_list_str[1:]):
                    #print("{}-{}".format(a, b))
                    bond_list.append("{}-{}".format(a, b))



        offset += N_b * bead.N
        #print("offset",offset)
        i += 1

    L = (system_mass/rho)**(1.0/3.0)


    snap = data.make_snapshot(N=N,
            box=data.boxdim(L=L),
            particle_types=particle_types,
            bond_types=bond_types)
    for _ in range(N):
        x = random.uniform(-L/2.0, L/2.0)
        y = random.uniform(-L/2.0, L/2.0)
        z = random.uniform(-L/2.0, L/2.0)
        point_list.append([x, y, z])

    snap.particles.typeid[:] = typeid_list
    snap.particles.mass[:] = mass_list
    snap.particles.position[:] = point_list
    for bond in bond_list:
        #print(bond)
        #print(bond.split("-")[0])
        #print(bond.split("-")[1])
        bond_a = int(bond.split("-")[0])
        bond_b = int(bond.split("-")[1])
        n_bonds = snap.bonds.N
        snap.bonds.resize(n_bonds + 1)
        snap.bonds.group[n_bonds] = [bond_a, bond_b]
        snap.bonds.typeid[n_bonds] = 0
    #print(bond_list)


    return snap



#print(snap.particles.position)
#print(snap.particles.mass)
#print(snap.particles.types)
#print(snap.particles.typeid[0])
#print(snap.particles.typeid[-1])
#print(snap.bonds.N)
#n_bonds = snap.bonds.N
#snap.bonds.resize(n_bonds + 1)
#print(snap.bonds.group)
#snap.bonds.group[n_bonds] = [0, 1]
#print(snap.bonds.typeid)
#print(snap.bonds.group)
#print(snap.bonds.types)
#n_bonds = snap.bonds.N
#snap.bonds.resize(n_bonds +1)
#snap.bonds.group[n_bonds] = [2, 3]
#print(snap.bonds.group)
#snap.bonds.typeid[0] = 0
#print(snap.bonds.N)

if __name__ == "__main__":
    from hoomd import md
    from hoomd import deprecated

    A = Bead()
    B = Bead(btype="B", mass = 1.0)
    C = PolyBead(btype="C", mass = 1.0, N = 10)
    snap = init_system({A : 100, B : 100, C : 5}, 0.1)


    init.read_snapshot(snap)
    deprecated.dump.xml(group = group.all(), filename = "start_init.hoomdxml", all=True)
    nl = md.nlist.cell()
    md.integrate.mode_standard(dt=2e-5)
    md.integrate.nve(group=group.all())
    dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=0.1, seed=0)
    dpd.pair_coeff.set('A', 'A', A=1.0, gamma = 1.0, r_cut = 1.0)
    dpd.pair_coeff.set('B', 'B', A=1.0, gamma = 1.0, r_cut = 1.0)
    dpd.pair_coeff.set('C', 'C', A=1.0, gamma = 1.0, r_cut = 1.0)

    dpd.pair_coeff.set('A', 'B', A=10.0, gamma = 1.0, r_cut = 1.0)
    dpd.pair_coeff.set('A', 'C', A=10.0, gamma = 1.0, r_cut = 1.0)
    dpd.pair_coeff.set('B', 'C', A=10.0, gamma = 1.0, r_cut = 1.0)

    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set('C-C', k=100.0, r0=1.0)
    analyze.log(filename="test_out.log", quantities=["pair_dpd_energy","volume","momentum","potential_energy","kinetic_energy","temperature","pressure", "bond_harmonic_energy"], period=1e4, header_prefix='#', overwrite=True)
    dump.dcd(filename="test_traj.dcd", period=5e4, overwrite=True)
    run(5e7)
    deprecated.dump.xml(group = group.all(), filename = "end_init.hoomdxml", all=True)
