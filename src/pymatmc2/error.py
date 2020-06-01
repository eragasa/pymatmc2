class Pymatmc2Error(Exception):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        try:
            self.msg = args[0]
        except IndexError as e:
            self.msg = ""

    def __str__(self):
        return self.msg

    def explain(self):
        return "{}: {}".format(
            self.__class__.__name__,
            self.msg)

class Pymatmc2ConcentrationMatrixError(Pymatmc2Error): pass
