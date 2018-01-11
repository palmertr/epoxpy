from enum import Enum


class Stages(Enum):
    MIXING = 1
    CURING = 2


class Integrators(Enum):
    NVE = 1
    NPT = 2
    LANGEVIN = 3

class NeighbourList(Enum):
    CELL = 1
    TREE = 2
