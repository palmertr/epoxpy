from epoxpy.tests.base_test import BaseTest
import pytest


class TestLegacyBonding(BaseTest):
    """
    Test class for testing simulation result for ABCTypeEpoxySimulation with legacy bonding enabled with baseline
    simulation result.
    """
    @pytest.mark.long
    @pytest.mark.frued_bonding
    def test_epoxy_sim_freud_bonding(self, datadir, tmpdir):
        """
        Checks if positions of particles are close to baseline particle positions. Here the baseline is legacy bonding
        where initial structure is not shrunk to desired density. Because the density is very less, the legacy bonding
        and freud bonding has the same behaviour and hence the trajectories are the same. But at higher density, due to
        the fact that freud neighbourlist returns neighbours sorted by distance and legacy neighbourlist returns
        neighbours sorted by particle id, the trajectory differs.
        :param datadir: 
        :param tmpdir: 
        :return: 
        """
        import epoxpy.abc_type_epoxy_simulation as es
        import epoxpy.job as jb
        import epoxpy.temperature_profile_builder as tpb
        import random
        import os
        import numpy as np
        import hoomd.data

        random.seed(1020)
        print('\n# Test: test_epoxy_sim_legacy_bonding')
        expected_gsd_file = os.path.join(datadir, 'legacy_bonding.gsd')
        print('expected gsd file path:{}'.format(expected_gsd_file))

        mix_time = 3e4
        mix_kt = 2.0
        time_scale = 100
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(250 * time_scale, 1.0 * temp_scale)

        out_dir = str(tmpdir)
        initial_structure_path = os.path.join(datadir, 'no_shrink_init.hoomdxml')
        myEpoxySim = es.ABCTypeEpoxySimulation('freud_bonding', mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, bond=True, n_mul=1.0,
                                               output_dir=out_dir, shrink=False,
                                               ext_init_struct_path=initial_structure_path)

        mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
        mySingleJobForEpoxy.execute()

        total_time = type_A_md_temp_profile.get_total_sim_time()
        gsd_write_period = myEpoxySim.dcd_write
        total_frames = int(round(total_time/gsd_write_period))
        print('total_frames:{}'.format(total_frames))

        current_gsd = tmpdir.join('freud_bonding', 'data.gsd')
        gsd_path = str(current_gsd)

        snapshot = hoomd.data.gsd_snapshot(gsd_path, total_frames-1)
        expected_snapshot = hoomd.data.gsd_snapshot(expected_gsd_file, total_frames-1)
        assert snapshot.particles.N == expected_snapshot.particles.N
        expected_pos = expected_snapshot.particles.position
        current_pos = snapshot.particles.position
        assert np.allclose(expected_pos, current_pos)
        expected_bonds = expected_snapshot.bonds.N
        current_bonds = snapshot.bonds.N
        print('test_epoxy_sim_freud_bonding_count. legacy bonds:{}, freud bonds:{}'.format(expected_bonds,
                                                                                           current_bonds))
        assert current_bonds == expected_bonds

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
        import epoxpy.abc_type_epoxy_simulation as es
        import epoxpy.job as jb
        import epoxpy.temperature_profile_builder as tpb
        import epoxpy.bonding as bondClass
        import random
        import os
        import hoomd.data
        import numpy as np

        random.seed(1020)
        expected_gsd_file = os.path.join(datadir, 'shrunk_freud_bonding.gsd')
        print('expected gsd file path:{}'.format(expected_gsd_file))

        mix_time = 3e4
        mix_kt = 2.0
        time_scale = 100
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(250 * time_scale, 1.0 * temp_scale)

        out_dir = str(tmpdir)
        sim_name = 'shrunk_freud_bonding'
        initial_structure_path = os.path.join(datadir, 'shrunk_init.hoomdxml')
        myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, bond=True, n_mul=1.0, shrink=True,
                                               output_dir=out_dir, ext_init_struct_path=initial_structure_path)

        mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
        mySingleJobForEpoxy.execute()

        total_time = type_A_md_temp_profile.get_total_sim_time()
        gsd_write_period = myEpoxySim.dcd_write
        total_frames = int(round(total_time/gsd_write_period))
        print('total_frames:{}'.format(total_frames))

        current_gsd = tmpdir.join(sim_name, 'data.gsd')
        gsd_path = str(current_gsd)

        snapshot = hoomd.data.gsd_snapshot(gsd_path, total_frames-1)
        expected_snapshot = hoomd.data.gsd_snapshot(expected_gsd_file, total_frames-1)
        assert snapshot.particles.N == expected_snapshot.particles.N
        expected_bonds = expected_snapshot.bonds.N
        current_bonds = snapshot.bonds.N
        print('test_epoxy_sim_freud_shrunk_regression. expected:{}, current:{}'.format(expected_bonds, current_bonds))
        assert current_bonds == expected_bonds
        expected_pos = expected_snapshot.particles.position
        current_pos = snapshot.particles.position
        assert np.allclose(current_pos, expected_pos)

        idxs, counts = np.unique(snapshot.bonds.group, return_counts=True)
        for idx, index in enumerate(idxs):
            p_typeid = snapshot.particles.typeid[idx]
            p_type = snapshot.particles.types[p_typeid]
            if p_type == 'A':
                assert (counts[index] <= bondClass.FreudBonding.MAX_A_BONDS)
            elif p_type == 'B':
                assert (counts[index] <= bondClass.FreudBonding.MAX_B_BONDS)
