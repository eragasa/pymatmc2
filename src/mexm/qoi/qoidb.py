import copy,yaml
from copy import deepcopy
from collections import OrderedDict
from pypospack.io.filesystem import OrderedDictYAMLLoader

class QoiConfiguration(object):
    """ Qoi Database

        Attributes:
            filename(str): file to read/write to yaml file
    """

    def __init__(self, qoi_dict=None):
        assert any([
            isinstance(qoi_dict, dict),
            qoi_dict is None
        ])

        self.filename = 'pypospack.qoi.yaml'
        self.qois = {}

        if qoi_dict is not None:
            self.configure_from_dict(qoidb=qoidb_OrderedDict)

    def initialize_from_dict(self, qoi_dict):
        o = QoiConfiguration()
        o.configure_from_dict(qoi_dict=qoi_dict)
        return o

    @property
    def qoi_names(self):
        return [k for k in self.qois.keys()]

    def configure_from_dict(self, qoi_dict):
        for qoi_k, qoi_info in qoi_dict.items():
            kwarg_names = ['qoi_name', 'qoi_type', 'structures', 'target', 'options']
            kwargs = {}
            for k in kwarg_names:
                try:
                    kwargs[k] = qoi_info[k]
                except KeyError as e:
                    if k == 'options':
                        pass
                    else:
                        raise

            self.add_qoi(**kwargs)

    def add_qoi(self,
            qoi_name,
            qoi_type,
            structures,
            target,
            options=None):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """
        assert isinstance(qoi_name,str)
        assert isinstance(qoi_type,str)
        assert isinstance(qoi_type, dict)
        assert isinstance(target, float)

        #<--------- create the dictionary entry for this qoi
        self.qois[qoi_name] = {
            'qoi_name':qoi_name,
            'qoi_type':qoi_type,
            'structures':deepcopy(structures),
            'target':target,
            'options':options
        }

    def read(self,filename=None):
        """ read qoi configuration from yaml file
        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """
        assert isinstance(filename,str)

        # set the attribute if not none
        if filename is not None:
            self.filename = filename

        try:
            with open(self.filename,'r') as f:
                self.qois = yaml.load(f, OrderedDictYAMLLoader)
        except:
            raise

        # <------------------
        self.qoi_names = [k for k in self.qois]

    def write(self,filename=None):
        """ write qoi configuration to yaml file

        Args:
            fname(str): file to write yaml from from.  If no argument is passed
               then use the filename attribute.  If the filename is set, then
               the filename attribute is also set.
        """

        # set the attribute if not none
        if filename is not None:
            self.filename = filename

        # marshall attributes into a dictionary
        _qoidb = copy.deepcopy(self.qois)

        # dump to yaml file
        with open(self.filename,'w') as f:
            yaml.dump(_qoidb,f, default_flow_style=False)

    def to_string(self):
        str_out = 80*'-'+"\n"
        str_out += '{:^80}\n'.format('QUANTITIES OF INTEREST')
        str_out += 80*'-'+"\n"
        for k,v in self.qois.items():
            str_out += '{}\n'.format(k)

        return str_out
#------------------------------------------------------------------------------
