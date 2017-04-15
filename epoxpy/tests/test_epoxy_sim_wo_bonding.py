from epoxpy.tests.base_test import BaseTest
import pytest


class TestWOBonding(BaseTest):
    """
    Test class for testing simulation result for ABCTypeEpoxySimulation with baseline simulation result.
    Checks if positions of particles are close to baseline particle positions.
    """
    @pytest.mark.long
    def test_epoxy_sim_wo_bonding(self, datadir, tmpdir):
        import epoxpy.abc_type_epoxy_simulation as es
        import epoxpy.job as jb
        import epoxpy.temperature_profile_builder as tpb
        import random
        import os
        import numpy as np
        import hoomd.data

        random.seed(1020)
        print('\n# Test: test_epoxy_sim_wo_bonding')
        expected_gsd_file = os.path.join(datadir, 'wo_bonding.gsd')
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
        exclude_mixing_in_output = False
        initial_structure_path = os.path.join(datadir, 'no_shrink_init.hoomdxml')
        myEpoxySim = es.ABCTypeEpoxySimulation('wo_bonding', mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, output_dir=out_dir, n_mul=1.0,
                                               exclude_mixing_in_output=exclude_mixing_in_output, shrink=False,
                                               ext_init_struct_path=initial_structure_path)

        mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
        mySingleJobForEpoxy.execute()

        total_time = type_A_md_temp_profile.get_total_sim_time()
        md_time = total_time-mix_time
        gsd_write_period = myEpoxySim.dcd_write
        total_frames = int(round(total_time/gsd_write_period))
        total_md_frames = int(round(md_time/gsd_write_period))
        print('total_frames:{}, md frames:{}'.format(total_frames, total_md_frames))

        print('current directory: {}'.format(os.getcwd()))
        print('tmp dir: {}'.format(tmpdir))
        print('datadir: {}'.format(datadir))
        print('expected gsd path:{}'.format(expected_gsd_file))
        #current_gsd = os.path.join(tmpdir,'wo_bonding/data.gsd')
        current_gsd = tmpdir.join('wo_bonding', 'data.gsd')

        gsd_path = str(current_gsd)
        if exclude_mixing_in_output is False:
            total_md_frames = total_frames
        snapshot = hoomd.data.gsd_snapshot(filename=gsd_path, frame=(total_md_frames-1))
        expected_snapshot = hoomd.data.gsd_snapshot(expected_gsd_file, total_frames-1)
        assert snapshot.particles.N == expected_snapshot.particles.N
        expected_pos = expected_snapshot.particles.position
        current_pos = snapshot.particles.position
        assert np.allclose(expected_pos, current_pos)

    @pytest.mark.long
    def test_epoxy_sim_wo_bonding_exclude_mixing_in_output(self, datadir, tmpdir):
        """
        Tests whether the exclude mixing in output breaks any existing functionality.
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
        print('\n# Test: test_epoxy_sim_wo_bonding')
        expected_gsd_file = os.path.join(datadir, 'exclude_mixing.gsd')
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
        exclude_mixing_in_output = True
        initial_structure_path = os.path.join(datadir, 'no_shrink_init.hoomdxml')
        myEpoxySim = es.ABCTypeEpoxySimulation('wo_bonding', mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, output_dir=out_dir, n_mul=1.0,
                                               exclude_mixing_in_output=exclude_mixing_in_output, shrink=False,
                                               ext_init_struct_path=initial_structure_path)

        mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
        mySingleJobForEpoxy.execute()

        total_time = type_A_md_temp_profile.get_total_sim_time()
        md_time = total_time-mix_time
        gsd_write_period = myEpoxySim.dcd_write
        total_frames = int(round(total_time/gsd_write_period))
        total_md_frames = int(round(md_time/gsd_write_period))
        print('total_frames:{}, md frames:{}'.format(total_frames, total_md_frames))

        print('current directory: {}'.format(os.getcwd()))
        print('tmp dir: {}'.format(tmpdir))
        print('datadir: {}'.format(datadir))
        print('expected gsd path:{}'.format(expected_gsd_file))
        #current_gsd = os.path.join(tmpdir,'wo_bonding/data.gsd')
        current_gsd = tmpdir.join('wo_bonding', 'data.gsd')

        gsd_path = str(current_gsd)
        if exclude_mixing_in_output is False:
            total_md_frames = total_frames
        snapshot = hoomd.data.gsd_snapshot(filename=gsd_path, frame=(total_md_frames-1))
        expected_snapshot = hoomd.data.gsd_snapshot(expected_gsd_file, frame=(total_md_frames-1))
        assert snapshot.particles.N == expected_snapshot.particles.N
        expected_pos = expected_snapshot.particles.position
        current_pos = snapshot.particles.position
        assert np.allclose(expected_pos, current_pos)
