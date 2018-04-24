import random
import math
import mbuild as mb
import numpy as np
from epoxpy.lib import A, B


class AB(mb.Compound):
    """A coarse grained dimer of one A and one B particle to allow
        hoomd to equilibrate system before beginning crosslinking routine"""
    
    mass = 2.0 
    
    def __init__(self, rotate_random=True):
        super(AB, self).__init__()
        num_particles = 2
                
        self.add(A(), label='A')
        self.add(B(), label='B')
        
        if rotate_random is True:
            mb.Compound.rotate(self['A'], random.uniform(0, 2*math.pi), around=np.asarray([1, 0, 0]))
            mb.Compound.rotate(self['A'], random.uniform(0, 2*math.pi), around=np.asarray([0, 1, 0]))
            mb.Compound.rotate(self['A'], random.uniform(0, 2*math.pi), around=np.asarray([0, 0, 1]))
        
        
        mb.force_overlap(move_this=self['B'],
                         from_positions=self['B']['down'],
                         to_positions=self['A']['up'])
