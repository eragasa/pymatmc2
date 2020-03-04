import os
from multicell import MultiCell
from pymatmc2 import Pymatmc2Configuration
if __name__ == "__main__":
    configuration = Pymatmc2Configuration()
    configuration.read(path=os.path.join('resource','pymatmc2.config'))

    simulation_type = 'vasp'
    simulation_paths = {
        'fcc1':os.path.join('simulations','fcc1_000_T400_P0'),
        'fcc2':os.path.join('simulations','fcc2_000_T400_P0')
    }

    obj = MultiCell.initialize_from_pymatmc2_configuration(
        configuration=configuration
    )

    obj.cells = None
    for simulation_name, simulation_path in simulation_paths.items():
        obj.add_simulation(
            simulation_type=simulation_type,
            simulation_name=simulation_name,
            simulation_path=simulation_path   
        )

    print(obj.cell_names)
    print(obj.cell_molar_fraction)
    print(obj.symbols)
    print(obj.total_molar_fraction)
    print(obj.phase_molar_fraction)