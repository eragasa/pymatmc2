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

    old_obj = MultiCell.initialize_from_pymatmc2_configuration(
        configuration=configuration
    )

    old_obj.cells = None
    for simulation_name, simulation_path in simulation_paths.items():
        
        old_obj.add_simulation(
            simulation_type=simulation_type,
            simulation_name=simulation_name,
            simulation_path=simulation_path   
        )

    print(old_obj.cell_names)
    print(old_obj.cell_molar_fraction)
    print(old_obj.symbols)
    print(old_obj.total_molar_fraction)
    print(old_obj.phase_molar_fraction)
    print(
        old_obj.simulations['fcc1'].poscar.path,
        [
            old_obj.simulations['fcc1'].poscar.get_number_of_atoms(s) for s in old_obj.symbols
        ],
        old_obj.simulations['fcc1'].outcar.total_energy
    )
    print(
        old_obj.simulations['fcc2'].poscar.path,
        [
            old_obj.simulations['fcc2'].poscar.get_number_of_atoms(s) for s in old_obj.symbols
        ],
        old_obj.simulations['fcc2'].outcar.total_energy
    )

    new_obj = MultiCell.initialize_from_pymatmc2_configuration(
        configuration=configuration
    )
    simulation_paths = {
        'fcc1':os.path.join('simulations','fcc1_001_T400_P0'),
        'fcc2':os.path.join('simulations','fcc2_001_T400_P0')
    }
    new_obj.cells = None
    for simulation_name, simulation_path in simulation_paths.items():
        new_obj.add_simulation(
            simulation_type=simulation_type,
            simulation_name=simulation_name,
            simulation_path=simulation_path   
        )
    print(new_obj.cell_names)
    print(new_obj.cell_molar_fraction)
    print(new_obj.symbols)
    print(new_obj.total_molar_fraction)
    print(new_obj.phase_molar_fraction)
    print(
        new_obj.simulations['fcc1'].poscar.path,
        [
            new_obj.simulations['fcc1'].poscar.get_number_of_atoms(s) for s in new_obj.symbols
        ],
        [
            new_obj.cells['fcc1'].get_number_of_atoms(s) for s in new_obj.symbols
        ],

        new_obj.simulations['fcc1'].outcar.total_energy
    )
    print(
        new_obj.simulations['fcc2'].poscar.path,
        [
            new_obj.simulations['fcc2'].poscar.get_number_of_atoms(s) for s in new_obj.symbols
        ],
        new_obj.simulations['fcc2'].outcar.total_energy
    )