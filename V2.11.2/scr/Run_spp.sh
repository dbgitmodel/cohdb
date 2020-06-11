#!/bin/sh

cd $PBS_O_WORKDIR

source /space/hpc-apps/profile.d/obs.sh
source /space/hpc-apps/profile.d/bira.sh
source /space/hpc-apps/profile.d/meteo.sh
module load mpt/2.12
module load 18/hdf-netcdf

time mpiexec_mpt  ./coherens

