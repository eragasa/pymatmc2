#! /bin/bash
#PBS -A PAA0028
#PBS -l walltime=3:59:00
#PBS -l nodes=8:ppn=28
#PBS -N fcc1_000_T400_P0 
#PBS -j oe

module load intel/19.0.5
module load intelmpi/2019.3

cd $PBS_O_WORKDIR

mpiexec $VASP_STD_BIN

exit

