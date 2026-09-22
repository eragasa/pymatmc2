from mexm.simulation.base import Simulation
from mexm.simulation.base import AtomicSimulation

from mexm.simulation.lammps import LammpsSimulation
from mexm.simulation.gulp import GulpSimulation
from mexm.simulation.phonts import PhontsSimulation
from mexm.simulation.vasp import VaspSimulation

from mexm.simulation.min_all import StructuralMinimization
from mexm.simulation.lmps_min_all import LammpsStructuralMinimization

class GulpStructuralMinimization(GulpSimulation, StructuralMinimization):
    simulation_type = 'gulp_min_all'

class VaspStructuralMinimization(VaspSimulation, StructuralMinimization):
    simulation_type = 'vasp_min_all'
    
class StaticCalculation(Simulation): pass

from mexm.simulation.lmps_min_none import LammpsStaticCalculation
class GulpStaticCalculation(GulpSimulation, StaticCalculation):
    simulation_type = 'gulp_min_none'
    
class VaspStaticCalculation(VaspSimulation, StaticCalculation):
    simulation_type = 'vasp_min_none'

class PositionMinimization(Simulation):
    pass

from mexm.simulation.lmps_min_pos import LammpsPositionMinimization

class GulpPositionMinimization(GulpSimulation, PositionMinimization):
    simulation_type = 'gulp_min_positions'

class VaspPositionMinimization(VaspSimulation, PositionMinimization):
    simulation_type = 'vasp_min_positions'
    
class PhononCalculation(Simulation):
    pass
class GulpPhononCalculation(GulpSimulation, PhononCalculation):
    simulation_type = 'gulp_phonons'
class VaspPhononCalculation(VaspSimulation, PhononCalculation):
    simulation_type = 'vasp_phonons'

class NptSimulation(Simulation): pass
from mexm.simulation.lmps_npt import LammpsNptSimulation

class NvtSimulation(Simulation): pass
class LammpsNvtSimulation(LammpsSimulation, NvtSimulation):
    simulation_type = 'lammps_nvt'

def all_subclasses(cls):
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses() for s in all_subclasses(c)]
    )

class SimulationFactory():
    factories = dict(
        **{k.simulation_type: k for k in LammpsSimulation.__subclasses__()},
        **{k.simulation_type: k for k in GulpSimulation.__subclasses__()},
        **{k.simulation_type: k for k in PhontsSimulation.__subclasses__()},
        **{k.simulation_type: k for k in VaspSimulation.__subclasses__()}
    )

    @staticmethod
    def create_Simulation(simulation_type: str) -> Simulation:
        return SimulationFactory.factories[simulation_type]()
#if False:


from mexm.simulation.manager import SimulationManager
from mexm.simulation.manager import SerialSimulationManager
from mexm.simulation.manager import MpiSimulationManager
