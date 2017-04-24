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

    return status


def run_epoxy_sim(sim_name, mix_time, mix_kt, temp_prof, bond, n_mul, shrink, legacy_bonding, ext_init_struct_path,
                  exclude_mixing_in_output, log_curing, curing_log_period, log_write, dcd_write, job, dt, density,
                  bond_period):

    fig_path = os.path.join(job.workspace(), 'temperature_profile.png')
    temp_temperature_profile = tpb.LinearTemperatureProfileBuilder(0)
    temp_temperature_profile.set_raw(temp_prof)
    temp_prof = temp_temperature_profile
    print('tempearture profile:{}'.format(temp_prof))
    fig = temp_prof.get_figure()
    fig.savefig(fig_path)
    in_path = os.path.join(job.workspace(), 'script_bckp.py')
    shutil.copy(__file__, in_path)

    myEpoxySim = es.ABCTypeEpoxySimulation(sim_name, mix_time=mix_time, mix_kt=mix_kt, temp_prof=temp_prof, bond=bond,
                                           n_mul=n_mul, shrink=shrink, legacy_bonding=legacy_bonding,
                                           ext_init_struct_path=ext_init_struct_path,
                                           exclude_mixing_in_output=exclude_mixing_in_output, log_curing=log_curing,
                                           curing_log_period=curing_log_period, log_write=log_write,
                                           dcd_write=dcd_write, output_dir=job.workspace(), dt=dt, density=density,
                                           bond_period=bond_period)

    mySingleJobForEpoxy = jb.SingleJob(myEpoxySim)
    mySingleJobForEpoxy.execute()

    job.document['bond_percent'] = myEpoxySim.get_curing_percentage()
    log_path = os.path.join(job.workspace(), 'curing.log')
    np.savetxt(log_path, myEpoxySim.curing_log)

    curing_log = zip(*myEpoxySim.curing_log)
    fig = plt.figure()
    plt.xlabel('Time steps')
    plt.ylabel('Cure percent')
    plt.margins(x=0.1, y=0.1)
    plt.plot(curing_log[0], curing_log[1])
    plt.plot(curing_log[0], curing_log[1], 'or')
    fig_path = os.path.join(job.workspace(), 'curing_curve.png')
    fig.savefig(fig_path)


def run_simulation(state_point):
    project = signac.init_project('ABCTypeEpoxy', 'data/')
    job = project.open_job(state_point)
    job.init()
    print('initialize', job)
    job_status = get_status(job)
    print('job status:{}'.format(job_status))
    if job_status == 'init':
        run_epoxy_sim(job=job, **job.statepoint())

kTs = [0.1]
mixing_temperature = 2.0
mixing_time = 3e4
time_scale = 1
for kT in kTs:
    flat_temp_profile = tpb.LinearTemperatureProfileBuilder(initial_temperature=mixing_temperature,
                                                            initial_time=mixing_time)
    flat_temp_profile.add_state_point(500 * time_scale, kT)
    sp = {'sim_name': 'epoxy_curing_flat_temperature_profile_{}kT'.format(kT),
          'mix_time': mixing_time,
          'mix_kt': mixing_temperature,
          'temp_prof': flat_temp_profile.get_raw(),
          'bond': True,
          'n_mul': 1.0,
          'shrink': True,
          'legacy_bonding': False,
          'ext_init_struct_path': None,
          'exclude_mixing_in_output': False,
          'log_curing': True,
          'curing_log_period': 1e5,
          'log_write': 1e5,
          'dcd_write': 1e5,
          'bond_period': 1e1,
          'dt': 1e-2,
          'density': 1.0}
    print(sp)
    run_simulation(sp)
