import mbuild as mb
from epoxpy.tests.base_test import BaseTest
import numpy as np


class TestBaseCompounds(BaseTest):

    def test_a(self, a):
        assert a.n_bonds == 0
        assert a.n_particles == 1

    def test_b(self, b):
        assert b.n_bonds == 0
        assert b.n_particles == 1

    def test_c(self, c):
        assert c.n_bonds == 0
        assert c.n_particles == 1

    def test_c10(self, c10):
        assert c10.n_bonds == 9
        assert c10.n_particles == 10

#    def test_c10_negative_spacing(self):
#        with pytest.raises(ValueError) as excinfo:
#            from epoxpy.lib import C10
#            C10(spacing=-1.0)
#        assert 'Molecule type C10 expects a positive non zero spacing.' in str(excinfo.value)

    def test_c10_position(self, c10, c10_ref_xyz):
        assert np.allclose(c10.xyz, c10_ref_xyz)

    def test_c10_new_position(self, c10_new_position, c10_new_ref_xyz):
        assert np.allclose(c10_new_position.xyz, c10_new_ref_xyz)

    def test_batch_add(self, a, b, c10):
        compound = mb.Compound()
        compound.add([a, b, c10])
        assert compound.n_particles == 1 + 1 + 10
        assert compound.n_bonds == 0 + 0 + 9

    def test_epoxy_blend_10a_20b_2c10(self, epoxy_a_10_b_20_c10_2_blend):
        assert epoxy_a_10_b_20_c10_2_blend.n_bonds == 18
        assert epoxy_a_10_b_20_c10_2_blend.n_particles == 50

    def test_epoxy_sim_wo_bonding(self):
        import epoxpy.epoxy_simulation as es
        import epoxpy.job as jb
        import epoxpy.temperature_profile_builder as tpb
        import random

        random.seed(1020)
        print('\n# Test1: Running the simulation in a single job')
        mix_time = 3e4
        mix_kt = 2.0
        time_scale = 1
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(240 * time_scale, 1.0 * temp_scale)

        myEpoxySim = es.EpoxySimulation('epoxy_test_mbuild', mix_time=mix_time, mix_kt=mix_kt,
                                        temp_prof=type_A_md_temp_profile, n_mul=1.0)

        mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
        mySingleJobForEpoxy.execute()
