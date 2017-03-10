from abc import ABCMeta, abstractmethod


class MdEngine(object):

    __metaclass__ = ABCMeta
    engine_name = 'None'

    def __init__(self, engine_name):
        self.engine_name = engine_name

    @abstractmethod
    def set_initial_structure(self, initialStructure):
        pass
