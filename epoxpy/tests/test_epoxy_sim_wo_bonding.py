from epoxpy.tests.base_test import BaseTest


class TestWOBonding(BaseTest):
    """
    Test class for testing simulation result for ABCTypeEpoxySimulation with baseline simulation result.
    Checks if positions of particles are close to baseline particle positions.
    """
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
        expected_gsd_file = os.path.join(datadir,'wo_bonding.gsd')
        print('expected gsd file path:{}'.format(expected_gsd_file))

        mix_time = 3e4
        mix_kt = 2.0
        time_scale = 100
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
        type_A_md_temp_profile.add_state_point(240 * time_scale, 1.0 * temp_scale)

        out_dir = str(tmpdir)
        myEpoxySim = es.ABCTypeEpoxySimulation('wo_bonding', mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, output_dir=out_dir, n_mul=1.0)

        mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
        mySingleJobForEpoxy.execute()

        total_time = type_A_md_temp_profile.get_total_sim_time()
        gsd_write_period = myEpoxySim.dcd_write
        total_frames = int(round(total_time/gsd_write_period))
        print('total_frames:{}'.format(total_frames))

        print('current directory: {}'.format(os.getcwd()))
        print('tmp dir: {}'.format(tmpdir))
        print('datadir: {}'.format(datadir))
        print('expected gsd path:{}'.format(expected_gsd_file))
        #current_gsd = os.path.join(tmpdir,'wo_bonding/data.gsd')
        current_gsd = tmpdir.join('wo_bonding','data.gsd')

        gsd_path = str(current_gsd)
        snapshot = hoomd.data.gsd_snapshot(filename=gsd_path, frame=(total_frames-1))
        expected_snapshot = hoomd.data.gsd_snapshot(expected_gsd_file, total_frames-1)
        assert snapshot.particles.N == expected_snapshot.particles.N
        expected_pos = expected_snapshot.particles.position
        current_pos = snapshot.particles.position
        assert np.allclose(expected_pos, current_pos)

