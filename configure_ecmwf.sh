#!/bin/bash

#COMPILER=intel
COMPILER=gcc

set -ex

if [ $COMPILER = intel ]
then
  module load prgenv/intel
  module load intel/2021.4.0
  module load intel-mpi/2021.4.0
  module load intel-mkl/19.0.5
  module load netcdf4/4.7.4
  module load hdf5/1.10.6
else
    module swap gcc/14.2.0
fi

# Compile options

# Default (optimized) settings
#CXXFLAGS="-Wall -g -O3 -march=native -std=c++11 -DADEPT_FAST_EXPONENTIAL -fopenmp"
CXXFLAGS="-Wall -g -O3 -march=native -std=c++11 -fopenmp"

# Optimized but with checking
#CXXFLAGS="-Wall -g -O2 -march=native -std=c++11 -DADEPT_BOUNDS_CHECKING -DADEPT_INIT_REAL_SNAN"

# Debug settings

#CXXFLAGS="-Wall -g -O0 -march=native -std=c++11 -DADEPT_BOUNDS_CHECKING -DADEPT_INIT_REAL_SNAN -fopenmp"
#CXXFLAGS="-Wall -g -O0 -march=native -std=c++11 -DADEPT_BOUNDS_CHECKING -fopenmp"

# Location of Adept automatic differentiation library
#ADEPT_VER=adept-2.1.3-intel
ADEPT_VER=adept-2.1.3-gcc14.2.0
ADEPT_DIR=/home/parr/apps/$ADEPT_VER
ADEPT_FLAGS="--with-adept=$ADEPT_DIR"

# Location of NetCDF-4 library
#module load netcdf4
#NETCDF_FLAGS="--with-netcdf=$NETCDF4_DIR"

#LDFLAGS=-Wl,-rpath,/usr/local/apps/szip/2.1/LP64/lib64

# Set install location
INSTALL_DIR=/home/parr/apps/tcrad-0.3.1-gcc14.2.0
#INSTALL_DIR=/home/parr/apps/tcrad-0.3-intel

# Call configure script
./configure --prefix "$INSTALL_DIR" "CXXFLAGS=$CXXFLAGS" $ADEPT_FLAGS $@
