#import matplotlib.pyplot as plt
import sys
import epoxpy.epoxy_simulation as es
import epoxpy.temperature_profile_builder as tpb


jobIndex = sys.argv[1]
jobs = []

# This simulation should run a total of 700 time steps because the default dt of the HOOMD engine is 1e-2
mix_time = 3e4
md_time = 4e4
mix_kt = 2.0
time_scale = 1
temp_scale = 1


type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_A_md_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
type_A_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
type_A_md_temp_profile.add_state_point(240 * time_scale, 1.0 * temp_scale)
fig = type_A_md_temp_profile.get_figure()
#plt.show()
fig.savefig("type_A_md_temp_profile.png")
myEpoxySim = es.EpoxySimulation('epoxy_test1', mix_time=mix_time, mix_kt=mix_kt, md_time=md_time,
                                temp_prof=type_A_md_temp_profile, n_mul=1.0)

type_B_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_B_md_temp_profile.add_state_point(30 * time_scale, 3.0 * temp_scale)
type_B_md_temp_profile.add_state_point(30 * time_scale, 4.5 * temp_scale)
type_B_md_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
type_B_md_temp_profile.add_state_point(240 * time_scale, 1.0 * temp_scale)
fig = type_B_md_temp_profile.get_figure()
#plt.show()
fig.savefig("type_B_md_temp_profile.png")
myEpoxySim2 = es.EpoxySimulation('epoxy_test2', mix_time=mix_time, mix_kt=mix_kt, md_time=md_time,
                                temp_prof=type_B_md_temp_profile, n_mul=1.0)


jobs.append(myEpoxySim)
jobs.append(myEpoxySim2)

jobs[jobIndex].execute()
