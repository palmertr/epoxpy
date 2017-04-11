from abc import ABCMeta, abstractmethod


class Job(object):

    __metaclass__ = ABCMeta
    sims = []

    def __init__(self):
        print('New job created')

    @abstractmethod
    def execute(self):
        pass


class BatchJob(Job):

    def __init__(self, simulation_objects):
        self.sims = simulation_objects
        print('New batch job created to run:')
        for sim in self.sims:
            print(sim.get_sim_name())

    def execute(self):
        for sim in self.sims:
            sim.execute()


class SingleJob(BatchJob):

    def __init__(self, simulation):
        del self.sims[:]
        self.sims.append(simulation)
        print('New single job created to run {}'.format(self.sims[0].get_sim_name()))
