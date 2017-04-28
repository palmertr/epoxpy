from abc import ABCMeta, abstractmethod


class Simulation(object):
    """Common base class for all simulations."""
    __metaclass__ = ABCMeta

    simulation_name = 'Blank Simulation'
    engine_name = 'None'

    def __init__(self, engine_name):
        self.engine_name = engine_name

    @abstractmethod
    def get_sim_name(self):
        pass

    @abstractmethod
    def initialize(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def execute(self):
        pass
