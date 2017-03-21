import numpy as np

import mbuild as mb


class B(mb.Compound):
    """An B particle with two ports one facing up and one facing down."""
    def __init__(self):
        super(B, self).__init__()
        self.add(mb.Particle(name='B'))

        self.add(mb.Port(anchor=self[0]), 'up')
        self['up'].translate([0, 0.07, 0])

        self.add(mb.Port(anchor=self[0]), 'down')
        self['down'].translate([0, -0.07, 0])
