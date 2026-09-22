from collections import OrderedDict
from mexm.qoi import DefectFormationEnergy

class LammpsDefectFormationEnergy(DefectFormationEnergy):
    qoi_type = 'lmps_defect'
    qois_calculated = ['lammps_E_formation']
    is_base_class = False

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures,dict)
        assert 'ideal' in structures
        assert 'defect' in structures

        DefectFormationEnergy.__init__(self,
                                       qoi_name=_qoi_name,
                                       qoi_type=_qoi_type,
                                       structures=_structures)

    def determine_simulations(self):

        # 1. minimize the prototype bulk structure
        ideal_structure_name = self.structures['ideal']
        ideal_simulation_type = 'lmps_min_all'
        ideal_simulation_name = "{}.{}".format(
                ideal_structure_name,
                ideal_simulation_type)
        self.add_simulation(
                simulation_id='ideal',
                simulation_name=ideal_simulation_name,
                simulation_type=ideal_simulation_type,
                simulation_structure=ideal_structure_name)

        # 2. position minimization of the defect structre
        # '{}.a11_min_all' is used as the length of the a1 vector
        defect_structure_name = self.structures['defect']
        defect_simulation_type = 'lmps_min_pos'
        defect_simulation_name = '{}.{}'.format(
                defect_structure_name,
                defect_simulation_type)
        self.add_simulation(
                simulation_id='defect',
                simulation_name=defect_simulation_name,
                simulation_type=defect_simulation_type,
                simulation_structure=defect_structure_name,
                bulk_structure_name=defect_structure_name)

    def calculate_qois(self,task_results):
        _prefix = '{}.{}'.format(
                self.structures['defect'],
                self.qoi_type)
        s_name_defect = self.structures['defect']
        s_name_bulk   = self.structures['ideal']

        e_defect = task_results[
                "{}.lmps_min_pos.toten".format(s_name_defect)]
        e_bulk = task_results[
                "{}.lmps_min_all.toten".format(s_name_bulk)]
        n_atoms_defect = task_results[
                "{}.lmps_min_pos.natoms".format(s_name_defect)]
        n_atoms_bulk = task_results[
                "{}.lmps_min_all.natoms".format(s_name_bulk)]
        e_f = e_defect - n_atoms_defect/n_atoms_bulk*e_bulk

        self.qois = OrderedDict()
        self.qois['{}.E_formation'.format(_prefix)]= e_f
