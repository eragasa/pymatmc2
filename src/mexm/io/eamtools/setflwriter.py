from mexm.io.eamtools import BaseSetflFile

class SetflWriter(BaseSetflFile):
    N_VALUES_PER_LINE_RHO = 5
    N_VALUES_PER_LINE_R = 5
    SETFL_NUM_FORMAT = "{:+24.16E}"
    SETFL_INT_FORMAT = "{:5d}"
    PAIR_KEY_FORMAT = "{}.{}"

    def __init__(self, path=None):
        BaseSetflFile.__init__(self, path=path)
        if self.path is not None:
            self.write(path=path)

    def write(self, path=None):
        if path is not None:
            assert isinstance(path, str)
            self.path = path
        raise NotImplementedError()
