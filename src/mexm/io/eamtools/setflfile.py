from mexm.io.eamtools import SetflReader 
from mexm.io.eamtools import SetflWriter

class SetflFile():
    """ class for the setfl

    """
    def __init__(self):
        self.src_file = SetflReader()
        self.dst_file = SetflWriter()

    @property
    def src_path(self):
        """:str path to read the setflfile"""
        return self.src_file.path

    @src_path.setter
    def src_path(self, path):
        self.src_file.path = path

    @property
    def dst_path(self):
        """:str path to write the setflfile"""
        return self.dst_file.path

    def read(self, path=None):
        """ read the setfl file

        Args:
            path (str): path to read the setflfile
        """

        # handle method arguments
        if path is not None:
            assert isinstance(path, str)
            self.src_path = path

        # read the file using the SetflReader
        self.src_file.read(path=self.src_path)

    def write(self, path=None):
        """ write the setfl file

        Args:
            path (str): path to write the setflfile
        """
        if path is not None:
            assert isinstance(path, None)
            self.dst_path = path
        self.dst_path.write(path=self.dst_path)
