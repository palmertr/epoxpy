import epoxpy.abc_type_epoxy_simulation as es
import epoxpy.job as jb
import epoxpy.temperature_profile_builder as tpb
import os
import numpy as np
import shutil
import matplotlib
import signac
matplotlib.use('webagg')
import matplotlib.pyplot as plt


def get_status(job):
    status = 'init'
    if job.isfile('final.hoomdxml') and job.isfile('out.log'):
        status = 'job-computed'
    elif 'temp_prof' in job.sp:
        status = 'temperature-profile-created'

    return status


def run_epoxy_sim(sim_name, mix_time, mix_kt, temp_prof, bond, n_mul, shrink, legacy_bonding, ext_init_struct_path,
                  exclude_mixing_in_output, log_curing, curing_log_period, log_write, dcd_write, job, dt, density,
                  bond_period, activation_energy, sec_bond_weight,use_dybond_plugin):
    fig_path = os.path.join(job.workspace(), 'temperature_profile.png')
    temp_temperature_profile = tpb.LinearTemperatureProfileBuilder(0)
    temp_temperature_profile.set_raw(temp_prof)
    temp_prof = temp_temperature_profile
    print('tempearture profile:{}'.format(temp_prof))
    #fig = temp_prof.get_figure()
    #fig.savefig(fig_path)
    in_path = os.path.join(job.workspace(), 'script_bckp.py')
    # shutil.copy(__file__, in_path)

    myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt, temp_prof=temp_prof, bond=bond,
                                           n_mul=n_mul, shrink=shrink, legacy_bonding=legacy_bonding,
                                           ext_init_struct_path=ext_init_struct_path,
                                           exclude_mixing_in_output=exclude_mixing_in_output, log_curing=log_curing,
                                           curing_log_period=curing_log_period, log_write=log_write,
                                           dcd_write=dcd_write, output_dir=job.workspace(), dt=dt, density=density,
                                           bond_period=bond_period, activation_energy=activation_energy,
                                           sec_bond_weight=sec_bond_weight,
                                           use_dybond_plugin=use_dybond_plugin)

    mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
    mySingleJobForEpoxy.execute()

    job.document['bond_percent'] = myEpoxySim.get_curing_percentage()
    log_path = os.path.join(job.workspace(), 'curing.log')
    np.savetxt(log_path, myEpoxySim.curing_log)
    bond_rank_log_path = os.path.join(job.workspace(), 'bond_rank.log')
    #print(myEpoxySim.bond_rank_log)
    np.savetxt(bond_rank_log_path,myEpoxySim.bond_rank_log)
    curing_log = list(zip(*myEpoxySim.curing_log))
    if len(curing_log) > 0:
        #fig = plt.figure()
        plt.xlabel('Time steps')
        plt.ylabel('Cure percent')
        plt.margins(x=0.1, y=0.1)
        plt.plot(curing_log[0], curing_log[1])
        plt.plot(curing_log[0], curing_log[1], 'or')
        fig_path = os.path.join(job.workspace(), 'curing_curve.png')
        #fig.savefig(fig_path)


def init_job(state_point):
    project = signac.init_project('ABCTypeEpoxy', 'data/')
    job = project.open_job(state_point)
    job.init()
    job_status = get_status(job)
    print('job status:{}'.format(job_status))
    if job_status == 'init':
        fig_path = os.path.join(job.workspace(), 'temperature_profile.png')
        temp_temperature_profile = tpb.LinearTemperatureProfileBuilder(0)
        temp_temperature_profile.set_raw(job.sp.temp_prof)
        temp_prof = temp_temperature_profile
        print('tempearture profile:{}'.format(temp_prof))
        #fig = temp_prof.get_figure()
        #fig.savefig(fig_path)

    print('initialize', job)
    return job


def run_simulation(state_point, Force=False):
    project = signac.init_project('ABCTypeEpoxy', 'data/')
    job = project.open_job(state_point)
    job.init()
    print('initialize', job)
    job_status = get_status(job)
    print('job status:{}'.format(job_status))
    if job_status == 'temperature-profile-created' or Force:
        run_epoxy_sim(job=job, **job.statepoint())


long_simulation = True

if long_simulation:
    time_scale = 10000
    n_mul = 1000.0
    curing_log_period = 1e5
else:
    time_scale = 10
    n_mul = 10.0
    curing_log_period = 1

kTs = [1.0]
mixing_temperature = 20.0
mixing_time = 3e4
jobs = []

for kT in kTs:
    flat_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mixing_temperature,
                                                            initial_time=mixing_time)
    flat_temp_profile.add_state_point(1 * time_scale, kT)
    flat_temp_profile.add_state_point(499 * time_scale, kT)

    sp = {'sim_name': 'epoxy_curing_flat_temperature_profile_{}kT'.format(kT),
          'mix_time': mixing_time,
          'mix_kt': mixing_temperature,
          'temp_prof': flat_temp_profile.get_raw(),
          'bond': True,
          'n_mul': n_mul,
          'shrink': True,
          'use_dybond_plugin':False,
          'legacy_bonding': True,
          'ext_init_struct_path': None,
          'exclude_mixing_in_output': False,
          'log_curing': True,
          'curing_log_period': curing_log_period,
          'log_write': 1,
          'dcd_write': 1,
          'bond_period': 1e1,
          'dt': 1e-2,
          'density': 1.0,
          'activation_energy': 1.0,
          'sec_bond_weight': 1.0}
    job = init_job(sp)
    jobs.append(job)

for job in jobs:
    run_simulation(job.statepoint(), Force=True)
