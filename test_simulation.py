import epoxpy.abc_type_epoxy_simulation as es
import epoxpy.job as jb
import epoxpy.temperature_profile_builder as tpb

mix_time = 3e4
mix_kt = 2.0
time_scale = 1
temp_scale = 1
cure_kt = 2.0
nmul = 1.0
log_period = 1e5
dump_period = 1e5
curing_log_period = 1e5

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
fig.savefig("type_A_temp_profile.png")

shrink = True
legacy_bonding = False
bonding = True
ext_init = False
exclude_mixing_in_output = False

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

myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                       temp_prof=type_A_md_temp_profile, bond=bonding, n_mul=nmul,
                                       shrink=shrink, legacy_bonding=legacy_bonding,
                                       ext_init_struct_path=initial_structure_path,
                                       exclude_mixing_in_output=exclude_mixing_in_output, log_curing=True,
                                       curing_log_period=curing_log_period,
                                       log_write=log_period,
                                       dcd_write=dump_period)

mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
mySingleJobForEpoxy.execute()

curing_log = zip(*myEpoxySim.curing_log)

snapshot = myEpoxySim.system.take_snapshot(bonds=True)
n_bonds = len(snapshot.bonds.group)
possible_bonds = ((10*4)+(20*2)+(2*9)) * nmul
bond_percent = (n_bonds/possible_bonds)*100

print('possible bonds:{}, made bonds:{}, percent cure:{}'.format(possible_bonds, n_bonds, bond_percent))

import matplotlib
matplotlib.use('webagg')
import matplotlib.pyplot as plt
import os
import numpy as np
import shutil

fig = plt.figure()
plt.xlabel('Time steps')
plt.ylabel('Cure percent')
plt.margins(x=0.1, y=0.1)
plt.plot(curing_log[0], curing_log[1])
plt.plot(curing_log[0], curing_log[1], 'or')
fig_path = os.path.join(sim_name, 'curing_curve_{}_kT.png'.format(cure_kt))
log_path = os.path.join(sim_name, 'curing.txt')
fig.savefig(fig_path)
np.savetxt(log_path, myEpoxySim.curing_log)

in_path = os.path.join(sim_name, 'script_bckp.py')
shutil.copy(__file__, in_path)
