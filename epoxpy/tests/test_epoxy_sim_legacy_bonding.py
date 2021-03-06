from epoxpy.tests.base_test import BaseTest
import pytest


class TestLegacyBonding(BaseTest):
    """
    Test class for testing simulation result for ABCTypeEpoxySimulation with legacy bonding enabled with baseline
    simulation result.
    Checks if positions of particles are close to baseline particle positions.
    """
    @pytest.mark.long
    def test_epoxy_sim_legacy_bonding(self, datadir, tmpdir):
        import epoxpy.abc_type_epoxy_dpd_simulation as es
        import epoxpy.temperature_profile_builder as tpb
        import random
        import os
        import numpy as np
        import gsd.hoomd

        random.seed(1020)


        shrink = False
        legacy_bonding = True
        bonding = True
        exclude_mixing_in_output = False
        mix_time = 3e4
        mix_kt = 2.0
        time_scale = 100
        cure_kt = 2.0
        nmul = 2.0
        log_period = 1e5
        dump_period = 1e2
        curing_log_period = 1e1

        type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
        type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)

        sim_name = 'legacy_bonding'
        out_dir = str(tmpdir)
        out_dir = os.path.join(out_dir, sim_name)
        myEpoxySim = es.ABCTypeEpoxyDPDSimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                               temp_prof=type_A_md_temp_profile, bond=bonding, n_mul=nmul,
                                               shrink=shrink, legacy_bonding=legacy_bonding,
                                               exclude_mixing_in_output=exclude_mixing_in_output, log_curing=False,
                                               curing_log_period=curing_log_period,
                                               log_write=log_period,
                                               dcd_write=dump_period,
                                               output_dir=out_dir)

        myEpoxySim.execute()

        current_gsd = tmpdir.join(sim_name, 'data.gsd')
        gsd_path = str(current_gsd)
        print('reading gsd: ', gsd_path)
        f = gsd.fl.GSDFile(gsd_path, 'rb')
        t = gsd.hoomd.HOOMDTrajectory(f)
        snapshot = t[-1]
        assert snapshot.particles.N == 100
        current_bonds = snapshot.bonds.N
        print('test_epoxy_sim_legacy_bonding_count. legacy bonds:{}'.format(current_bonds))
        assert current_bonds > 30#Just checking if some bonds are being made
