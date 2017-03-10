from factories import SimulationEngineFactory
from abc import ABCMeta, abstractmethod


class Simulation(object):
    """Common base class for all simulations."""
    __metaclass__ = ABCMeta

    simulation_name = 'Blank Simulation'
    engine_name = 'None'

    def __init__(self, engine_name):
        self.engine = SimulationEngineFactory.get_engine(engine_name)

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
