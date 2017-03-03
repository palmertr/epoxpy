from hoomd_engine import HOOMDEngine


class SimulationEngineFactory(object):

    @classmethod
    def get_engine(cls, engine_name):
        engine = 'None'
        if engine_name == 'HOOMD':
            engine = HOOMDEngine()

        if engine == 'None':
            raise ValueError('Specified simulation engine not made by this factory!')

        print('{} created in SimulationEngineFactory'.format(engine))
        return engine


class JobSchedulerFactory(object):

    @classmethod
    def slurm():
        return 'slurm'

    @classmethod
    def get_scheduler(cls, scheduler_name):
        scheduler = 'none'
        if scheduler_name == 'HOOMD':
            scheduler = 'hoomd_engine'
        print('{} created in JobSchedulerFactory'.format(scheduler))
        return scheduler
