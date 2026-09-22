import os
import shutil

from mexm.simulation import AtomicSimulation
from mexm.io.vasp import Incar
from mexm.io.vasp import Poscar
from mexm.io.vasp import Potcar
from mexm.io.vasp import Kpoints
from mexm.io.vasp import Outcar
from mexm.io.vasp import Oszicar
from mexm.io.vasp import Contcar

class VaspSimulation(AtomicSimulation):
    """

    Attributes:
        incar(Incar)
        poscar(Poscar)
        kpoints(Kpoints)
        outcar(Outcar): initialized as None, object initialized on read()
        contcar(Contcar): initailized as None, object initialized on read()
        oszicar(Ozicar): initialized as None, object initialized on read()
    """

    def __init__(self):
        self.src_path = None
        self.dst_path = None
        self.incar = Incar()
        self.poscar = Poscar()
        self.potcar = Potcar()
        self.kpoints = Kpoints()
        self.outcar = None
        self.contcar = None
        self.oszicar = None

    @property
    def structure(self) -> Poscar:
        if isinstance(self.contcar, Contcar):
            return self.contcar
        else:
            return self.poscar

    @property
    def total_energy(self):
        return self.outcar.total_energy

    def write(self, simulation_path):
        self.dst_path = simulation_path

        poscar_path = os.path.join(simulation_path, 'POSCAR')
        potcar_path = os.path.join(simulation_path, 'POTCAR')
        incar_path = os.path.join(simulation_path, 'INCAR')
        kpoints_path = os.path.join(simulation_path, 'KPOINTS')

        self.poscar.write(path=poscar_path)
        self.potcar.symbols = self.poscar.symbols
        self.potcar.write(path=potcar_path)
        self.incar.write(path=incar_path)
        self.kpoints.write(path=kpoints_path)

    def read(self, simulation_path):
        self.src_path = simulation_path

        # these files are the input files
        poscar_path = os.path.join(simulation_path, 'POSCAR')
        potcar_path = os.path.join(simulation_path, 'POTCAR')
        incar_path = os.path.join(simulation_path, 'INCAR')
        kpoints_path = os.path.join(simulation_path, 'KPOINTS')

        self.poscar.read(path=poscar_path)
        self.incar.read(path=incar_path)
        self.kpoints.read(path=kpoints_path)
        
        try:
            self.potcar.read(path=potcar_path)
        except FileNotFoundError:
            pass

        # these files are the output files
        oszicar_path = os.path.join(simulation_path, 'OSZICAR')
        outcar_path = os.path.join(simulation_path, 'OUTCAR')
        contcar_path = os.path.join(simulation_path, 'CONTCAR')

        if os.path.isfile(oszicar_path):
            self.oszicar = Oszicar()
            self.oszicar.read(path=oszicar_path)

        if os.path.isfile(outcar_path):
            self.outcar = Outcar()
            self.outcar.read(path=outcar_path)

        if os.path.isfile(contcar_path):
            self.contcar = Contcar()
            self.contcar.read(path=contcar_path)
    
    def update_hpc_configuration(self, n_cores, n_nodes):
        pass

    def archive(
        self, 
        dst_path: str,
        files_to_archive=[
            'POSCAR', 'INCAR', 'KPOINTS', 
            'OSZICAR', 'OUTCAR', 'CONTCAR', 'vasprun.xml'
        ]
    ):
        """
        Arguments:
            dst_path (str): directory which the VASP files will be created
            files_to_archive (List[str]): which VASP files to archive, by default
                this method will try to archive the POSCAR, INCAR, KPOINTS, 
                OSZICAR, OUTCAR, CONTCAR, and vasprun.xml 
        """

        self.dst_path = dst_path

        # delete destination path if it exists
        if os.path.isdir(self.dst_path):
            shutil.rmtree(self.dst_path)
        os.mkdir(self.dst_path)

        # archive the files        
        for f in files_to_archive:
            try:
                src_file_path = os.path.join(self.src_path, f)
                dst_file_path = os.path.join(self.dst_path, f)
                shutil.copy(src = src_file_path, dst = dst_file_path)
            except FileNotFoundError:
                pass
