#!/bin/bash

echo ""
echo "Note: This install script assumes: "
echo "1. Availibility of C++ compilers with c++11 support"
echo "2. Availability of MPI compilers"
echo ""

# Cleanup any previous installation

## ./uninstall.sh


# Grab the specific stable branch of LAMMPS compaitbility has been tested for


if [[ ! -d build/lammps_stable_29Oct2020 ]] ; then
    echo "Cloning lammps_stable code in the build directory"    
    mkdir -p build/lammps_stable_29Oct2020
    git clone --depth 1 --branch stable_29Oct2020 https://github.com/lammps/lammps.git build/lammps_stable_29Oct2020
else
    echo "Using existing lammps_stable code in the build directory"
fi   


# Copy ChIMES files to correct locations

## Use a pre-built ChIMES JIT MODEL
if [[ -z "$1" ]] ; then
    echo "Usage: install_jit.sh <model>"
    exit 1
else
    MODEL=$1
fi

## Double-check that required files are available.
if [[ -f ../../chimesFF/src/chimesJIT.h ]] ; then
    echo "Found chimesJIT.h"
else
    echo "Did not find chimesJIT.h"
    exit 1
fi
if [[ -f ../../chimesFF/src/chimesJIT.cpp ]] ; then
    echo "Found chimesJIT.cpp"
else
    echo "Did not find chimesJIT.cpp"    
    exit 1
fi
if [[ -f ../../chimesFF/src/chimesJITev.$MODEL.h ]] ; then
    echo "Found ../../chimesFF/src/chimesJITev.$MODEL.h"
else
    echo "Did not find ../../chimesFF/src/chimesJITev.$MODEL.h"    
    exit 1
fi
if [[ -f ../../chimesFF/src/chimesJITev.$MODEL.cpp ]] ; then
    echo "Found ../../chimesFF/src/chimesJITev.$MODEL.cpp"
else
    echo "Did not find ../../chimesFF/src/chimesJITev.$MODEL.cpp"    
    exit 1
fi
if [[ -f ../../chimesFF/src/chimesJITcluster.$MODEL.cpp ]] ; then
    echo "Found ../../chimesFF/src/chimesJITcluster.$MODEL.cpp"
else
    echo "Did not find ../../chimesFF/src/chimesJITcluster.$MODEL.cpp"    
    exit 1
fi

echo "BUILDING A JIT VERSION OF LAMMPS FOR $MODEL"

cp ../../chimesFF/src/chimesJIT.{h,cpp}	build/lammps_stable_29Oct2020/src/MANYBODY/
cp ../../chimesFF/src/chimesJITev.$MODEL.h	build/lammps_stable_29Oct2020/src/MANYBODY/chimesJITev.h
cp ../../chimesFF/src/chimesJITev.$MODEL.cpp	build/lammps_stable_29Oct2020/src/MANYBODY/chimesJITev.cpp
cp ../../chimesFF/src/chimesJITcluster.$MODEL.cpp	build/lammps_stable_29Oct2020/src/MANYBODY/chimesJITcluster.cpp
cp src/pair_chimes_jit.{h,cpp} 		build/lammps_stable_29Oct2020/src/MANYBODY/
cp etc/pair.{h,cpp} 			build/lammps_stable_29Oct2020/src
cp etc/Makefile.mpi_chimes 		build/lammps_stable_29Oct2020/src/MAKE

# Compile

cd build/lammps_stable_29Oct2020/src
make yes-manybody
make mpi_chimes
cd -

# Finish
if [[ ! -d exe ]] ; then
    mkdir exe
fi    
mv build/lammps_stable_29Oct2020/src/lmp_mpi_chimes exe/lmp_mpi_chimes.$MODEL

loc=`pwd`
echo ""
echo "Compilation complete. "
echo "Generated the following LAMMPS executable with ChimesJIT support:"
echo "${loc}/exe/lmp_mpi_chimes.$MODEL"
echo "See ${loc}/tests for usage examples"
echo ""



