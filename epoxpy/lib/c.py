import mbuild as mb


class C(mb.Compound):
    """An coarse grained particle called C(arbitrarily) with two ports one facing up and one facing down."""
    def __init__(self):
        super(C, self).__init__()

        self.add(mb.Particle(name='C'))

        self.add(mb.Port(anchor=self[0]), 'up')
        mb.translate(self['up'], [0, 0.07, 0])

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.translate(self['down'], [0, -0.07, 0])
