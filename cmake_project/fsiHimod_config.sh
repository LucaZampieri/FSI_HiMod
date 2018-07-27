#!/bin/bash
#
# Example configuration file for LifeV
#
# Andrea Bartezzaghi, 14-Nov-2017
#


echo "Entered fsiHimod_config.sh "

# compilers to use
export C_COMPILER=gcc
export CXX_COMPILER=g++
export FORTRAN_COMPILER=gfortran

export MPI_C_COMPILER=mpicc
export MPI_CXX_COMPILER=mpic++
export MPI_FORTRAN_COMPILER=mpif90
export MPI_EXEC=mpiexec

# Src directory (src?)
FSIHIMOD=src

# additional postfix to add to the build directory
POSTFIX=

# build type, can either be "Release" or "Debug"
BUILD_TYPE=Release
#BUILD_TYPE=Debug

# module selection
# for each module, uncomment and set to ON to force its inclusion, uncomment and
# set to OFF to disable it, leave commented to use the default setting

ENABLE_ETA=ON

# set this to pass additional parameters to cmake
PARAMS=

echo "Exiting fsiHimod_config.sh "
