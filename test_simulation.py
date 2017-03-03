from epoxy_simulation import EpoxySimulation
from job import SingleJob
from temperature_profile_builder import LinearTemperatureProfileBuilder

print('\n# Test1: Running the simulation in a single job')
# This simulation should run a total of 700 time steps because the default dt of the HOOMD engine is 1e-2
mix_time = 3e4
md_time = 4e4
mix_kt = 1.0
time_scale = 1
temp_scale = 1
type_A_temp_profile = LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_A_temp_profile.add_state_point(60 * time_scale, 4.5 * temp_scale)
type_A_temp_profile.add_state_point(190 * time_scale, 4.5 * temp_scale)
type_A_temp_profile.add_state_point(240 * time_scale, 1.0 * temp_scale)

myEpoxySim = EpoxySimulation('epoxy_test1', mix_time=mix_time, mix_kt=mix_kt, md_time=md_time,
                             temp_prof=type_A_temp_profile, n_mul=1.0)
mySingleJobForEpoxy = SingleJob(myEpoxySim)
mySingleJobForEpoxy.execute()
