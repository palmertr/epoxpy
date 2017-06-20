import epoxpy.resume_abc_type_epoxy_simulation as es
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
mix_time = 0
mix_kt = 2.0
time_scale = 100
temp_scale = 1
cure_kt = 2.0
nmul = 1.0
log_period = 1e5
dump_period = 1e2
curing_log_period = 1e1
cure_time = (500*time_scale)+mix_time
type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=cure_kt,
                                    initial_time=cure_time)

cool_down_temps = np.linspace(cure_kt,0.4,10)
print('cool_down_temps',cool_down_temps)
for temp in cool_down_temps:
    type_A_md_temp_profile.add_state_point(1, temp)
    type_A_md_temp_profile.add_state_point(100 * time_scale, temp)

fig = type_A_md_temp_profile.get_figure()
fig.savefig('myfigure.png')

shrink = False
legacy_bonding = True
bonding = False
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

#initial_structure_path=None
out_dir = str(tmpdir)
out_dir = os.path.join(out_dir, sim_name)
print('out_dir', out_dir)
myEpoxySim = es.ResumeABCTypeEpoxySimulation(sim_name,mix_kt=mix_kt,mix_time=mix_time,
                                       temp_prof=type_A_md_temp_profile, bond=bonding, n_mul=nmul,
                                       legacy_bonding=legacy_bonding,
                                       ext_init_struct_path=initial_structure_path,
                                       exclude_mixing_in_output=exclude_mixing_in_output, log_curing=False,
                                       curing_log_period=curing_log_period,
                                       log_write=log_period,
                                       dcd_write=dump_period,
                                       output_dir=out_dir,
                                       reset_random_after_initialize=True)


mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
mySingleJobForEpoxy.execute()

