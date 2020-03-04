#PBS -A PAA0028
#PBS -l walltime=4:00:00
#PBS -l nodes=4:ppn=40
#PBS -N fcc1_001_T400_P0
#PBS -e job.err
#PBS -o job.out
#PBS -S /bin/bash

cd $PBS_O_WORKDIR


module load intel/19.0.5
module load intelmpi/2019.3


mpiexec $VASP_STD_BIN > vasp.out
touch jobComplete
