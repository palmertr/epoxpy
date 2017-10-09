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
        import epoxpy.temperature_profile_builder as tpb
        import random
        import os
        import numpy as np
        import gsd.hoomd

        random.seed(1020)
        print('\n# Test: test_epoxy_sim_wo_bonding')

        mix_time = 3e4
        mix_kt = 2.0
        cure_kt = 2.0
        time_scale = 100
        temp_scale = 1
        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)

        sim_name = 'wo_bonding'
        out_dir = str(tmpdir)
        exclude_mixing_in_output = False
        out_dir = os.path.join(out_dir, sim_name)
        myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, output_dir=out_dir, n_mul=2.0,
                                               exclude_mixing_in_output=exclude_mixing_in_output, shrink=False)

        myEpoxySim.execute()

        total_time = type_A_md_temp_profile.get_total_sim_time()
        md_time = total_time-mix_time
        gsd_write_period = myEpoxySim.dcd_write
        total_frames = int(round(total_time/gsd_write_period))
        total_md_frames = int(round(md_time/gsd_write_period))
        print('total_frames:{}, md frames:{}'.format(total_frames, total_md_frames))

        print('current directory: {}'.format(os.getcwd()))
        print('tmp dir: {}'.format(tmpdir))
        print('datadir: {}'.format(datadir))

        current_gsd = tmpdir.join(sim_name, 'data.gsd')
        gsd_path = str(current_gsd)
        f = gsd.fl.GSDFile(gsd_path, 'rb')
        t = gsd.hoomd.HOOMDTrajectory(f)
        snapshot = t[-1]

        assert snapshot.particles.N == 100
