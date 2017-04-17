import epoxpy.abc_type_epoxy_simulation as es
import epoxpy.job as jb
import epoxpy.temperature_profile_builder as tpb
import random
print('\n# Test1: Running the simulation in a single job')
# This simulation should run a total of 700 time steps because the default dt of the HOOMD engine is 1e-2

random.seed(1020)
mix_time = 3e4
mix_kt = 2.0
time_scale = 100
temp_scale = 1
type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
type_A_md_temp_profile.add_state_point(250 * time_scale, 1.0 * temp_scale)

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
ext_init = True
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

myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                       temp_prof=type_A_md_temp_profile, bond=bonding, n_mul=1.0,
                                       shrink=shrink, legacy_bonding=legacy_bonding,
                                       ext_init_struct_path=initial_structure_path,
                                       exclude_mixing_in_output=exclude_mixing_in_output)

mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
mySingleJobForEpoxy.execute()
