import copy,importlib
from collections import OrderedDict

class QoiManager(object):

    def __init__(self, configuration=None, fullauto=True):
        assert any([
            configuration is None,
            isinstance(configuration, dict)
        ])

        self.db = {}
        self.qoi_map = {}
        self.simulations = {}

        self.configure_from_dict(configuration=configuration)

        if fullauto:
            self.configure()
            self.determine_simulations()

    @staticmethod
    def get_qoi_map():
        qoi_map = {}
        return qoi_map

    def add_qoi_object():
        module_name = self.qoi_map[qoi_type]['module']
        module_


    @property
    def qoi_names(self):
        return [k for k in self.db]

    def configure(self):
        if len(self.qoi_map) > 0:
            self.qoi_map = QoiManager.get_qoi_map()

    def configure_qoi_objects(self):
        for qoik, qoiv in self.db.items():
            assert qoi_ID

    @staticmethod
    def get_qoi_name(qoi_type, structures):
        """ get qoi name from standard QOI naming convention

        In order to coordinate the transfer of data between the QoiManager and
        the SimulationManager, it is necessary to standardize the naming
        conventions.  The user may use their own Qoi_IDs, but internally
        MEXM standardizes the naming convention by the `qoi_type` attribute
        combined with one or more structureID/structures associated with the
        simulation.

        Args:
            qoi_type (str)
            structures (str, list, dict)
        """
        assert isinstance(qoi_type, str)
        assert any([
            isinstance(structures, str),
            isinstance(structures, list),
            isinstance(structures, dict)
        ])

        qoi_name = self.qoi_map[qoi_type].get_qoi_name(structures)
        return qoi_name
