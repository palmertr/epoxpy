from abc import ABCMeta, abstractmethod
from hoomd import variant


class TemperatureProfileBuilder(object):
    """Common base class for all TemperatureProfiles."""
    __metaclass__ = ABCMeta

    def __init__(self):
        self.temperature_profile = []

    @abstractmethod
    def get_profile(self):
        """
            Abstract method for TemperatureProfileBuilder.get_profile()
            Returns a hoomd.variant object
        """
        pass


class LinearTemperatureProfileBuilder(TemperatureProfileBuilder):
    """Builds a Linear Temperature Profile."""
    def __init__(self, initial_temperature, initial_time=0):
        TemperatureProfileBuilder.__init__(self)
        self.initial_temperature = initial_temperature
        self.temperature_profile.append((initial_time, initial_temperature))

    def add_state_point(self, ramp_time, desired_temperature):
        """
            Add a state point for the linear temperature profile. Adds ramp time to previous time or initial_time
        """
        print('Adding state point to LinearTemperatureProfileBuilder')
        last_state_point = self.temperature_profile[-1]
        new_state_point = ((last_state_point[0]+ramp_time), desired_temperature)
        if new_state_point[0] > last_state_point[0]:
            self.temperature_profile.append(new_state_point)
        else:
            err_string = 'Inconsistent state point added. Previous time: {}, new time: {}'.format(last_state_point[0], new_state_point[0])
            raise ValueError(err_string)
        print('Added new state point to LinearTemperatureProfile: {}'.format(self.temperature_profile))

    def get_profile(self):
        """
            Returns a hoomd.variant object
        """
        return variant.linear_interp(points=self.temperature_profile)
