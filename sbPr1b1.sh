#!/bin/bash
#SBATCH -o Pr1b1.out
#SBATCH -N 1
#SBATCH -e slurm.err
#SBATCH -n 20
#SBATCH --mem-per-cpu=500MB
#SBATCH --time 144:00:00
module purge
module load OpenMPI
export kLengthSZ=107



mpifort  -O0 -fdefault-real-8 -fdefault-double-8 -fcheck=all -fcheck=bounds -fcheck=array-temps \
        -fcheck=pointer -ffpe-trap=invalid,zero,overflow \
        -finit-real=snan -finit-integer=-999 \
        -Wall -Wextra -fbacktrace -DKLENGTH_SIZE=$kLengthSZ array_dimensions.F90 kernel_interp.F90 integrationModules.F90 EDQNMstratifiedModules.F90 job_parameters.F90 getForcing.F90 data_loader_all.F90 timeRoutines.F90 mainSST.F90 -o mainSSTProgram.x



#mpifort  -O0  -fdefault-real-8 -fdefault-double-8 -fbacktrace kernel_interp.F90 -DKLENGTH_SIZE=$kLengthSZ array_dimensions.F90 integrationModules.F90 EDQNMstratifiedModules.F90 job_parameters.F90 getForcing.F90 data_loader_all.F90 timeRoutines.F90 mainSST.F90 -o mainSSTProgram.x

mpirun -n $SLURM_NTASKS ./mainSSTProgram.x
