import random
import numpy as np
#from freud import parallel
from freud import box, locality
from abc import ABCMeta, abstractmethod


class Bonding(object):
    __metaclass__ = ABCMeta

    MAX_A_BONDS = 4
    MAX_B_BONDS = 2

    def __init__(self, system, groups, log, cut_off_dist=1.0, activation_energy=0.1):
        self.system = system
        self.groups = groups
        self.log = log
        self.activation_energy = activation_energy
        #parallel.setNumThreads(4)
        self.rank_dict = {}
        self.get_rank = self.rank_dict.get
        self.cut_off_dist = cut_off_dist
        self.bond_rank_log = []

    def __call__(self, step):
        self.find_pair(step)

    @abstractmethod
    def find_pair(self, timestep):
        pass

    @staticmethod
    def make_bond(index_a, index_b, snapshot):
        n_bonds = snapshot.bonds.N
        snapshot.bonds.resize(n_bonds + 1)
        snapshot.bonds.group[n_bonds] = [index_a, index_b]
        # sets new bond to be A-B type
        snapshot.bonds.typeid[n_bonds] = 1  # we know A-B bond type's id is 1

    @staticmethod
    def get_bond_rank(index, snapshot):
        return np.count_nonzero(snapshot.bonds.group == index)

    @staticmethod
    def pbc_diff(p1, p2, axes):
        dr = p1 - p2
        for i, p in enumerate(dr):
            dr[i] = abs(dr[i])
            if dr[i] > axes[i]*0.5:
                dr[i] -= axes[i]
        return np.sqrt(dr[0]**2+dr[1]**2+dr[2]**2)

    @staticmethod
    def bond_test(kT, delta_e, bond_rank, sec_bond_weight=500):
        # No idea if this is thread safe, watch out for MPI gotchas
        # Should be able to to tune rate with delta_e and kT
        # delta_e = 1
        mb_stats = np.exp(-delta_e / kT)
        # Divides by bond rank to make it less probable, add one to prevent rank 0
        # issues
        weight = 1
        if bond_rank >= 1:
            weight = sec_bond_weight
        if mb_stats / float(weight) > random.random():
            return True
        else:
            return False


class LegacyBonding(Bonding):
    def find_neighbours_and_bond(self, axis, bond_from_idx, bond_from_type, bond_to_group, bond_to_max_rank,
                                 bond_to_type, snapshot, xyz0):
        made_bonds = False
        for p in bond_to_group:
            xyz1 = p.position
            r = self.pbc_diff(np.array(xyz0), np.array(xyz1), axis)
            if r < self.cut_off_dist:
                bond_to_idx = p.tag
                bond_to_rank = self.get_rank(bond_to_idx, -1)
                if bond_to_rank < 0:  # -1 is the default value returned when id is not in dict.
                    bond_to_rank = self.get_bond_rank(bond_to_idx, snapshot)
                    self.rank_dict[bond_to_idx] = 0

                if bond_to_rank < bond_to_max_rank:
                    if self.bond_test(self.log.query('temperature'), self.activation_energy, bond_to_rank):
                        # bond_test
                        self.make_bond(bond_from_idx, bond_to_idx, snapshot)
                        made_bonds = True
                        self.rank_dict[bond_from_idx] += 1
                        self.rank_dict[bond_to_idx] += 1
                        break

        return made_bonds

    def find_pair(self, timestep):
        # Until I can hack a way to get access to the neighborlist
        snapshot = self.system.take_snapshot(bonds=True)
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
            bond_to_group = self.groups[1]#.group_b
            bond_to_type = 'B'
            bond_to_max_rank = self.MAX_B_BONDS
            bond_from_max_rank = self.MAX_A_BONDS
            # typeB = "B"
            # MAX_RANK_B = self.MAX_B_BONDS
        elif bond_from_type == 'B':
            bond_to_group = self.groups[0]#.group_a
            bond_to_type = 'A'
            bond_to_max_rank = self.MAX_A_BONDS
            bond_from_max_rank = self.MAX_B_BONDS
            # typeB = "A"
            # MAX_RANK_B = self.MAX_A_BONDS
        # Check to see if it can make more bonds

        bond_from_rank = self.get_rank(bond_from_idx, -1)
        if bond_from_rank < 0:  # -1 is the default value returned when id is not in dict.
            bond_from_rank = self.get_bond_rank(bond_from_idx, snapshot)
            self.rank_dict[bond_from_idx] = 0

        if bond_from_rank < bond_from_max_rank:
            xyz0 = snapshot.particles.position[bond_from_idx]
            axis = [snapshot.box.Lx, snapshot.box.Ly, snapshot.box.Lz]

        #           start_time = time.time()
            made_bonds = self.find_neighbours_and_bond(axis, bond_from_idx, bond_from_type, bond_to_group,
                                                       bond_to_max_rank, bond_to_type, snapshot, xyz0)
                #           stop_time = time.time()
            # print("time to calc neighbours for 1 frame = {}".format(stop_time - start_time))
            if made_bonds is True:
                self.system.restore_snapshot(snapshot)


class FreudBonding(Bonding):
    def __init__(self, system, groups, log, activation_energy):
        Bonding.__init__(self, system=system, groups=groups, log=log, activation_energy=activation_energy)
        # create freud nearest neighbor object
        # set number of neighbors
        self.n_neigh = 6
        # create freud nearest neighbors object
        self.nn = locality.NearestNeighbors(rmax=self.cut_off_dist, n_neigh=self.n_neigh, strict_cut=True)
        snapshot = self.system.take_snapshot()
        self.fbox = box.Box(Lx=snapshot.box.Lx, Ly=snapshot.box.Ly, Lz=snapshot.box.Lz)

    def find_neighbours_and_bond(self, bond_from_idx, bond_from_type, bond_to_max_rank, bond_to_type, snapshot,
                                 time_step):
        made_bonds = False
        # compute nearest neighbors for 6 nearest neighbors
        self.nn.compute(self.fbox, snapshot.particles.position, snapshot.particles.position)
        # get the neighborlist
        n_list = self.nn.getNeighborList()
        n_idxs = n_list[bond_from_idx]
        for n_idx in n_idxs:
            if n_idx != self.nn.getUINTMAX():
                p_typeid = snapshot.particles.typeid[n_idx]
                p_type = snapshot.particles.types[p_typeid]
                if p_type == bond_to_type:
                    bond_to_rank = self.get_rank(n_idx, -1)
                    if bond_to_rank < 0:  # -1 is the default value returned when id is not in dict.
                        bond_to_rank = 0#self.get_bond_rank(n_idx, snapshot)
                        self.rank_dict[n_idx] = 0

                    if bond_to_rank < bond_to_max_rank:
                        if self.bond_test(self.log.query('temperature'), self.activation_energy, bond_to_rank):
                            self.make_bond(bond_from_idx, n_idx, snapshot)
                            made_bonds = True
                            self.rank_dict[bond_from_idx] += 1
                            self.rank_dict[n_idx] += 1
                            break
        return made_bonds

    def find_pair(self, timestep):
        # Until I can hack a way to get access to the neighborlist
        snapshot = self.system.take_snapshot(bonds=True)
        n_p = snapshot.particles.N
        bond_from_type = "C"
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
            bond_to_type = 'B'
            bond_to_max_rank = self.MAX_B_BONDS
            bond_from_max_rank = self.MAX_A_BONDS
            # typeB = "B"
            # MAX_RANK_B = self.MAX_B_BONDS
        elif bond_from_type == 'B':
            bond_to_type = 'A'
            bond_to_max_rank = self.MAX_A_BONDS
            bond_from_max_rank = self.MAX_B_BONDS

        bond_from_rank = self.get_rank(bond_from_idx, -1)
        if bond_from_rank < 0:  # -1 is the default value returned when id is not in dict.
            bond_from_rank = 0#self.get_bond_rank(bond_from_idx, snapshot)
            self.rank_dict[bond_from_idx] = 0

        if bond_from_rank < bond_from_max_rank:
            made_bonds = self.find_neighbours_and_bond(bond_from_idx, bond_from_type, bond_to_max_rank, bond_to_type,
                                                       snapshot,timestep)
            if made_bonds is True:
                self.system.restore_snapshot(snapshot)
