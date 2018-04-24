
# coding: utf-8

# In[3]:


import hoomd
import epoxpy.abc_NP_type_epoxy_lj_harmonic_simulation as es
import epoxpy.temperature_profile_builder as tpb
import epoxpy.bonding as bondClass
import random
import os
import gsd.hoomd
import numpy as np

hoomd.context.initialize("--mode=gpu")

random.seed(1020)

mix_time = 3e5
mix_kt = 2.0
cure_kt = 5.0
time_scale = 5000
temp_scale = 1

t_Final = 7e6
t_SetT = 0.75e6
curing_log_period = 1e5

activation_energy = 2.0
sec_bond_weight = 2.0 #twice the primary bond energy for secondary bond
stop_bonding_after_percent = 100.0

av_calibrationT = 439.3636#K
av_aij_A_B = 30.729
bonding_period = 10
percent_bonds_per_step = 0.05
kT =500/av_calibrationT
AB_bond_const=100.0
AB_bond_dist=1.0

type_A_md_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mix_kt, initial_time=mix_time)
type_A_md_temp_profile.add_state_point(500 * time_scale, cure_kt)

out_dir = str('.')
sim_name = 'ljNP'
out_dir = os.path.join(out_dir, sim_name)
myEpoxySim = es.ABCNPTypeEpoxyLJHarmonicSimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt,
                                       temp_prof=type_A_md_temp_profile,
                                       bond=True, n_mul=1.0, shrink=True, 
                                       num_a=5000, num_b=5000, num_spheres=5, radius_sphere=5.0, radius_a=0.5, radius_b=0.5,
                                       shrink_time=1e5,
                                       output_dir=out_dir,
                                       use_dybond_plugin=True,
                                       bond_period=bonding_period,
                                       activation_energy=activation_energy,
                                       sec_bond_weight=sec_bond_weight,
                                       stop_after_percent=stop_bonding_after_percent,
                                       percent_bonds_per_step=percent_bonds_per_step,
                                       AA_interaction=25.00,
                                       AB_interaction=av_aij_A_B,
                                       AB_bond_const=AB_bond_const,
                                       AB_bond_dist= AB_bond_dist)

myEpoxySim.execute()



# In[4]:


import matplotlib.pyplot
import numpy as np
data = np.loadtxt("./ljNP/out.log")
AB_conversion = data[:, -1]
print(AB_conversion)

