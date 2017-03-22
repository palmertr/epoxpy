import mbuild as mb
from epoxpy.tests.base_test import BaseTest


class TestBaseCompounds(BaseTest):

    def test_batch_add(self, a, b):
        compound = mb.Compound()
        compound.add([a, b])
        assert compound.n_particles == 1 + 1
        assert compound.n_bonds == 0 + 0
