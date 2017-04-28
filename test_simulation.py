import epoxpy.abc_type_epoxy_simulation as es
import epoxpy.job as jb
import epoxpy.temperature_profile_builder as tpb
import os
import numpy as np
import shutil
import matplotlib
import gsd
import gsd.fl
import gsd.hoomd
matplotlib.use('webagg')
import matplotlib.pyplot as plt
import random

random.seed(1020)
mix_time = 3e4
mix_kt = 2.0
time_scale = 100
temp_scale = 1
cure_kt = 2.0
nmul = 1.0
log_period = 1e5
dump_period = 1e2
curing_log_period = 1e1

type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)
#type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
#type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
#type_A_md_temp_profile.add_state_point(250 * time_scale, 1.0 * temp_scale)

# to show the plot in an interactive mode, uncomment the below. Otherwise, just look at the saved image of the
# temperature profile.
#import matplotlib
#matplotlib.use('qt4agg')
#import matplotlib.pyplot as plt
fig = type_A_md_temp_profile.get_figure()
#plt.show()


shrink = False
legacy_bonding = True
bonding = True
ext_init = False
exclude_mixing_in_output = True

if ext_init is True:
    if shrink is True:
        initial_structure_path = 'shrunk_init.hoomdxml'
    else:
        initial_structure_path = 'no_shrink_init.hoomdxml'
else:
    initial_structure_path = None

if bonding is False:
    sim_name = 'no_bonding'
else:
    if legacy_bonding is True:
        sim_name = 'legacy_bonding'
    else:
        sim_name = 'freud_bonding'

sim_name = '{}_ts_{}_nmul_{}'.format(sim_name, time_scale, nmul)
#sim_name = 'benchmark'
print('sim_name', sim_name)
fig_path = os.path.join(sim_name, 'type_A_temp_profile.png')

if not os.path.exists(sim_name):
    os.makedirs(sim_name)

#fig.savefig(fig_path)
in_path = os.path.join(sim_name, 'script_bckp.py')
shutil.copy(__file__, in_path)

tmpdir = '.'
datadir = './benchmark'
expected_gsd_file = os.path.join(datadir, 'data.gsd')
print('expected gsd file path:{}'.format(expected_gsd_file))

shrink = True
legacy_bonding = False
bonding = True
exclude_mixing_in_output = False
mix_time = 3e4
mix_kt = 2.0
time_scale = 100
cure_kt = 2.0
nmul = 1.0
log_period = 1e5
dump_period = 1e2
curing_log_period = 1e1

type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)


out_dir = str(tmpdir)
initial_structure_path = os.path.join(datadir, 'initial.hoomdxml')
#initial_structure_path=None
out_dir = os.path.join(out_dir, sim_name)
print('out_dir', out_dir)
myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                       temp_prof=type_A_md_temp_profile, bond=bonding, n_mul=nmul,
                                       shrink=shrink, legacy_bonding=legacy_bonding,
                                       ext_init_struct_path=initial_structure_path,
                                       exclude_mixing_in_output=exclude_mixing_in_output, log_curing=False,
                                       curing_log_period=curing_log_period,
                                       log_write=log_period,
                                       dcd_write=dump_period,
                                       output_dir=out_dir,
                                       reset_random_after_initialize=True)

mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
mySingleJobForEpoxy.execute()

total_time = type_A_md_temp_profile.get_total_sim_time()
gsd_write_period = myEpoxySim.dcd_write
total_frames = int(round(total_time/gsd_write_period))
print('total_frames:{}'.format(total_frames))

current_gsd = os.path.join(tmpdir, sim_name, 'data.gsd')
gsd_path = str(current_gsd)
print('reading gsd: ', gsd_path)
f = gsd.fl.GSDFile(gsd_path, 'rb')
t = gsd.hoomd.HOOMDTrajectory(f)
snapshot = t[-1]
f = gsd.fl.GSDFile(expected_gsd_file, 'rb')
t = gsd.hoomd.HOOMDTrajectory(f)
expected_snapshot = t[-1]
assert snapshot.particles.N == expected_snapshot.particles.N
expected_pos = expected_snapshot.particles.position
current_pos = snapshot.particles.position
assert np.allclose(expected_pos, current_pos)
expected_bonds = expected_snapshot.bonds.N
current_bonds = snapshot.bonds.N
print('test_epoxy_sim_legacy_bonding_count. benchmark:{}, current simulation:{}'.format(expected_bonds,
                                                                                    current_bonds))
assert current_bonds == expected_bonds

