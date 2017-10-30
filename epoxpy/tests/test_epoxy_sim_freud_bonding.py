from epoxpy.tests.base_test import BaseTest
import pytest


class TestLegacyBonding(BaseTest):
    """
    Test class for testing simulation result for ABCTypeEpoxySimulation with legacy bonding enabled with baseline
    simulation result.
    """
    @pytest.mark.long
    @pytest.mark.frued_bonding
    def test_epoxy_sim_freud_shrunk_regression(self, datadir, tmpdir):
        """
        Here we are doing regression testing for the new bonding routine that operates and the mbuild initial structure
        whose volume is shrunk to a density of 1.0
        :param datadir:
        :param tmpdir:
        :return:
        """
        import epoxpy.abc_type_epoxy_dpd_simulation as es
        import epoxpy.temperature_profile_builder as tpb
        import epoxpy.bonding as bondClass
        import random
        import os
        import gsd.hoomd
        import numpy as np

        random.seed(1020)

        mix_time = 3e4
        mix_kt = 2.0
        cure_kt = 2.0
        time_scale = 100
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)

        out_dir = str(tmpdir)
        sim_name = 'shrunk_freud_bonding'
        out_dir = os.path.join(out_dir, sim_name)
        myEpoxySim = es.ABCTypeEpoxyDPDSimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile,
                                               bond=True, n_mul=2.0, shrink=True,
                                               output_dir=out_dir,
                                               use_dybond_plugin=False)

        myEpoxySim.execute()

        current_gsd = tmpdir.join(sim_name, 'data.gsd')
        gsd_path = str(current_gsd)
        print('reading gsd: ', gsd_path)
        f = gsd.fl.GSDFile(gsd_path, 'rb')
        t = gsd.hoomd.HOOMDTrajectory(f)
        snapshot = t[-1]
        current_bonds = snapshot.bonds.N
        assert snapshot.particles.N == 100
        print('test_epoxy_sim_freud_shrunk_regression. current:{}'.format(current_bonds))
        assert current_bonds > 30 #Just checking if some bonds are being made

        idxs, counts = np.unique(snapshot.bonds.group, return_counts=True)
        for index,idx in enumerate(idxs):
            p_typeid = snapshot.particles.typeid[idx]
            p_type = snapshot.particles.types[p_typeid]
            if p_type == 'A':
                assert (counts[index] <= bondClass.FreudBonding.MAX_A_BONDS)
            elif p_type == 'B':
                assert (counts[index] <= bondClass.FreudBonding.MAX_B_BONDS)


    def test_total_possible_bonds(self, tmpdir, datadir):
        import epoxpy.abc_type_epoxy_dpd_simulation as es
        import epoxpy.temperature_profile_builder as tpb
        import os

        mix_time = 3e4
        mix_kt = 2.0
        cure_kt = 2.0
        time_scale = 100
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)

        out_dir = str(tmpdir)
        sim_name = 'shrunk_freud_bonding'
        out_dir = os.path.join(out_dir, sim_name)
        initial_structure_path = os.path.join(datadir, 'shrunk_init.hoomdxml')
        myEpoxySim = es.ABCTypeEpoxyDPDSimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                            temp_prof=type_A_md_temp_profile, bond=True, n_mul=1000.0, shrink=True,
                                            output_dir=out_dir, ext_init_struct_path=initial_structure_path,
                                            reset_random_after_initialize=True)
        possible_bonds = myEpoxySim.total_possible_bonds()
        expected_possible_bonds = 40000 #becuase n_mul is 1000.0 and 10,000 A particles and 20,000 B particles can make a maximum of 40,000 bonds
        assert(possible_bonds == expected_possible_bonds)
        curing_percent = myEpoxySim.get_curing_percentage()
        expected_curing_percent = 0
        assert(curing_percent == expected_curing_percent)
