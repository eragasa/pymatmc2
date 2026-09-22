import os, sys
from copy import deepcopy
# import bibtexparser

class ReferenceManager(object):
    """ Reference manager

    Args:
        reference_paths (list of str): list of paths to bibtex formatted files
    """

    def __init__(self, reference_paths=None):
        assert any([
            reference_paths is None,
            isinstance(reference_paths, str),
            isinstance(reference_paths, list)
        ])

        self.db = {} # dictionary of reference entries
        self.initialize_reference_paths(reference_paths=reference_paths)
        self.read_reference_paths()

    def initialize_reference_paths(self, reference_paths):
        self.reference_paths = [self.__get_base_reference_path()]

        if isinstance(reference_paths, str):
            reference_paths.append(str)
        elif isinstance(reference_paths, list):
            assert all([isinstance(k, str) for k in reference_paths])
            self.reference_paths += deecopy(reference_paths)

    def __get_base_reference_path(self):
        base_reference_dir = os.path.dirname(
            os.path.abspath(
                sys.modules[ReferenceManager.__module__].__file__
            )
        )
        base_reference_path = os.path.join(base_reference_dir, 'references.bib')
        return base_reference_path

    def read_reference_paths(self):
        self.db = {}
        for ref_path in self.reference_paths:
            with open(ref_path) as bibtex_file:
                db = bibtexparser.load(bibtex_file)
            for article in db.entries:
                self.db[article['ID']] = deepcopy(article)
