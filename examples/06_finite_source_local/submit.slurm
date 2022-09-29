#!/bin/bash --login

#SBATCH --job-name=JOB_NAME_XXX
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

#SBATCH --account=ACCOUNT_NAME_XXX
#SBATCH --partition=highmem
#SBATCH --qos=highmem

cd WORKING_DIRECTORY_ABSOLUTE_PATH

export OMP_NUM_THREADS=1

module load PrgEnv-gnu
module load cray-fftw/3.3.8.9
module load metis/5.1.0
module load cray-hdf5/1.12.0.3
module load cray-netcdf/4.7.4.3
module load cmake/3.21.3

srun --distribution=block:block --hint=nomultithread ./axisem3d > axisem.out

