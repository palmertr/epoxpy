import hoomd
import hoomd.md
import hoomd.data as hd
import hoomd.deprecated
import numpy.random as rd
import math
import numpy as np
import os
import sys

def callback1(timestep):
    #x = timestep+1
    snapshot = system.take_snapshot(bonds=True)
    N_p = snapshot.particles.N
    system.restore_snapshot(snapshot)

if __name__ == "__main__":
	callback = int(sys.argv[1])
	if callback==1:
		bond_freq = float(sys.argv[2])
		print('bonding frequency: ',bond_freq)

	kT=10
	cells_on_box_edge=24
	time_steps=1e5
	random_initial_structure=False, 
	init_temp=False
	tau=0.1
	log_period=1e3
	dt=0.0001
	job=None
	a = 1.72
	box_edge = a * cells_on_box_edge
	n_cells = cells_on_box_edge ** 3
	n_particles = 4 * n_cells  # 4 because we assume FCC


	log_path = 'out.log'
	gsd_path = 'data.gsd'
	print('log path:{}'.format(log_path))

	print('kT: {}'.format(kT))
	hoomd.context.initialize('')
	if random_initial_structure == True:
		system = hoomd.deprecated.init.create_random(N=n_particles, box=hd.boxdim(L=box_edge * (1.0)))
	else:
		system = hoomd.init.create_lattice(unitcell=hoomd.lattice.fcc(a=1.72), n=cells_on_box_edge)
	if init_temp is True:
		snapshot = system.take_snapshot()
		snapshot = initialize_snapshot_T(snapshot,kT)
		#hoomd.context.initialize('--mode=cpu')
		system.restore_snapshot(snapshot)
		#system = hoomd.init.read_snapshot(snapshot)
		#hoomd.run(1)

	snapshot1 = system.take_snapshot()
	v = snapshot1.particles.velocity
	KE = 1/2 * np.mean(v**2)
	print('KE:{}'.format(KE))
	#T = 2/3 * KE
	T = 2 * KE
	print('T:{}'.format(T))

	nl = hoomd.md.nlist.cell()
	all = hoomd.group.all()
	lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
	lj.pair_coeff.set('A', 'A', epsilon=1, sigma=1)
	hoomd.md.integrate.mode_standard(dt=0.0001)
	nvt = hoomd.md.integrate.nvt(group=all, kT=kT, tau=tau)
	quantities_to_log = ['potential_energy', 'kinetic_energy', 'temperature', 'pressure', 'ndof']

	hoomd.analyze.log(filename=log_path, quantities=quantities_to_log, period=log_period,
		      overwrite=True)
	hoomd.dump.gsd(gsd_path, period=log_period, group=all, overwrite=True)
	if callback==1:
		bond_callback = hoomd.analyze.callback(callback = callback1, period = bond_freq)
	hoomd.util.quiet_status()
	hoomd.run(time_steps)
	nvt.disable()

