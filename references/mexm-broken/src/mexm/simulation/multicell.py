from mexm.simulation import AtomicSimulation
from mexm.simulation import MexmSimulationException

class MultiCellSimulation(AtomicSimulation):
    is_base_class = True

    def __init__(self,
                 name,
                 simulation_path):
        AtomicSimulation.__init__(self,
                                  name = name,
                                  simulation_path = simulation_path)
        self.simulations = {}

    def add_simulation(self, 
                       simulation_name, 
                       simulation):
        """
        Args:
            simulation_name (str)
            simulation (AtomicSimulation
        """
        if not isinstance(simulation, AtomicSimulation):
            msg
            raise TypeError(msg)
        
        self.simulation[simulation_name] = {}
        self.simulation[simulation_name]['obj'] = simulation

    def create_simulation_path(self, path = None):
        AtomicSimulation.create_simulation_path(path=path)

        for k,v in self.simulations.simulation_name:
            

