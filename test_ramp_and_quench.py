#import matplotlib.pyplot as plt
import epoxpy.epoxy_simulation as es
import epoxpy.job as jb
import epoxpy.temperature_profile_builder as tpb

print('\n# Test1: Running the simulation in a single job')
# This simulation should run a total of 700 time steps because the default dt of the HOOMD engine is 1e-2
mix_time = 1e5
mix_kt = 2.0
time_scale = 1e2
temp_scale = 1
type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=1, initial_time=mix_time)
type_A_md_temp_profile.add_state_point(5 * time_scale, 1.0 * temp_scale)
type_A_md_temp_profile.add_state_point(1 * time_scale, 2.0 * temp_scale)
type_A_md_temp_profile.add_state_point(5 * time_scale, 2.0 * temp_scale)
type_A_md_temp_profile.add_state_point(1 * time_scale, 1.5 * temp_scale)
type_A_md_temp_profile.add_state_point(0 * time_scale, 0.5 * temp_scale)
type_A_md_temp_profile.add_state_point(5 * time_scale, 0.5 * temp_scale)

fig = type_A_md_temp_profile.get_figure()
#plt.show()
fig.savefig("type_A_temp_profile.png")

myEpoxySim = es.EpoxySimulation('epoxy_test_ramp_and_quench1', mix_time=mix_time, mix_kt=mix_kt,
                                temp_prof=type_A_md_temp_profile, n_mul=10.0)

mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
mySingleJobForEpoxy.execute()
