import hoomd
hoomd.context.initialize('')
system = hoomd.init.create_lattice(unitcell=hoomd.lattice.sq(a=1.2), n=16)
from bokeh.plotting import figure, output_file, show
uc = hoomd.lattice.unitcell(N=2,
a1=[2, 0, 0],
a2=[0, 2, 0],
a3=[0, 0, 1],
dimensions=2,
position=[[0,0,0], [0.6, 0.6, 0]],                                                                                                                              type_name=['A', 'B'])
# Get a snapshot from the unitcell and add in bond topology
snap = uc.get_snapshot()
snap.bonds.resize(1)
snap.bonds.group[0] = [0, 1]
snap.bonds.types = ['bondA']
 
# replicate the lattice and initialize hoomd
snap = system.take_snapshot()
p = figure(width=400, height=400)
p.circle(x=snap.particles.position[:,0], y=snap.particles.position[:,1], radius=0.5, radius_units="data")
show(p)
                   
