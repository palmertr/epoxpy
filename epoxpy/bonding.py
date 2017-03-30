import random
import numpy as np

class Bonding:
    CUT = 1.0
    MAX_A_BONDS = 4
    MAX_B_BONDS = 2

    def __init__(self, system, epoxy_sim, log):
        self.system = system
        self.epoxy_sim = epoxy_sim
        self.log = log

    def __call__(self, step):
        self.find_pair(step)

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
            bond_to_group = self.epoxy_sim.group_b
            bond_to_type = 'B'
            bond_to_max_rank = self.MAX_B_BONDS
            bond_from_max_rank = self.MAX_A_BONDS
            # typeB = "B"
            # MAX_RANK_B = self.MAX_B_BONDS
        elif bond_from_type == 'B':
            bond_to_group = self.epoxy_sim.group_a
            bond_to_type = 'A'
            bond_to_max_rank = self.MAX_A_BONDS
            bond_from_max_rank = self.MAX_B_BONDS
            # typeB = "A"
            # MAX_RANK_B = self.MAX_A_BONDS
        # Check to see if it can make more bonds

        bond_from_rank = self.get_bond_rank(bond_from_idx, snapshot)

        if bond_from_rank < bond_from_max_rank:
            xyz0 = snapshot.particles.position[bond_from_idx]
            axis = [snapshot.box.Lx, snapshot.box.Ly, snapshot.box.Lz]
            for p in bond_to_group:
                xyz1 = p.position
                r = self.pbc_diff(np.array(xyz0), np.array(xyz1), axis)
                if r < self.CUT:
                    bond_to_idx = p.tag
                    bond_to_rank = self.get_bond_rank(bond_to_idx, snapshot)
                    if bond_to_rank < bond_to_max_rank:
                        delta_e = 0.1
                        if self.bond_test(self.log.query('temperature'), delta_e, bond_to_rank):
                            # bond_test
                            self.make_bond(bond_from_idx, bond_to_idx, snapshot)
                            print("Found one, bonding {} ({}) to {} ({})".format(bond_from_type, bond_from_idx,
                                                                                 bond_to_type, bond_to_idx))
                            # print("Rank of A {} type of A {}".format(rank, typeA))
                            # print("Rank of B {} type of B {}".format(bond_rank, typeB))
            self.system.restore_snapshot(snapshot)
        else:
            print("Can't bond it anymore")

    @staticmethod
    def make_bond(index_a, index_b, snapshot):
        n_bonds = snapshot.bonds.N
        snapshot.bonds.resize(n_bonds + 1)
        snapshot.bonds.group[n_bonds] = [index_a, index_b]
        # sets new bond to be A-B type
        snapshot.bonds.typeid[n_bonds] = 1  # we know A-B bond type's id is 1

    @staticmethod
    def get_bond_rank(index, snapshot):
        rank = 0
        for bond in snapshot.bonds.group:
            if bond[0] == index or bond[1] == index:
                rank += 1
        return rank

    @staticmethod
    def pbc_diff(p1, p2, axes):
        dr = p1 - p2
        for i, p in enumerate(dr):
            dr[i] = abs(dr[i])
            if dr[i] > axes[i]*0.5:
                dr[i] -= axes[i]
        return np.sqrt(dr[0]**2+dr[1]**2+dr[2]**2)

    @staticmethod
    def bond_test(kT, delta_e, bond_rank):
        # No idea if this is thread safe, watch out for MPI gotchas
        # Should be able to to tune rate with delta_e and kT
        # delta_e = 1
        mb_stats = np.exp(-delta_e / kT)
        # Divides by bond rank to make it less probable, add one to prevent rank 0
        # issues
        weight = 1
        if bond_rank >= 1:
            weight = 500
        if mb_stats / float(weight) > random.random():
            return True
        else:
            return False
