"""
mexm exception classes
"""

class MexmException(Exception):
    def __init__(self, msg, original_exception=None):
        super(MexmException, self).__init__(msg)
        self.original_exception = original_exception

class BaseError(Exception):
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

# deprecated
class BaseException(BaseError): pass

#TODO: rneame *Exception names to *Error names

class MexmUnsupportedOperationSystemError(BaseError): pass

class EamSetflException(BaseException): pass

class MexmBadEamEosError(BaseException): pass

class BadPreprocessorTypeException(BaseException): pass

class PyposmatBadParameterError(BaseException): pass

class PyposmatSamplingTypeError(BaseException): pass

class BadParameterException(PyposmatBadParameterError): pass

class BadManifoldTypeException(BaseException): pass

class BadClusterTypeException(BaseException): pass

class BadNearestNeighborTypeException(BaseException): pass

class LammpsSimulationError(BaseException): pass

class PyposmatError(BaseException): pass

class PypospackTaskManagerError(BaseException): pass

class PypospackBadKdeBandwidthType(BaseException): pass

class PypospackUnknownDistributionType(BaseException): pass

class PyposmatUnknownDistributionType(BaseException): pass

class PyposmatUnknownQoiFilterType(BaseException): pass

class PypospackBadEamEosError(BaseException): pass
