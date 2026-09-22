from mpi4py import MPI
import sys

vasp_args = ['vasp']
comm = MPI.COMM_SELF.Spawn(sys.executable, args=vasp_args)
