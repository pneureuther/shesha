"""
COMPASS simulation package
Abstraction layer for initialization and execution of a COMPASS simulation
"""
__all__ = ['simulator', 'simulatorBrama', 'simulatorCACAO', 'bench', 'benchBrama']

from .simulator import Simulator
from .simulatorBrama import SimulatorBrama
from .simulatorCACAO import SimulatorCACAO
from .bench import Bench
from .benchBrama import BenchBrama
