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

cp ../../chimesFF/src/chimesFF.{h,cpp}	build/lammps_stable_29Oct2020/src/MANYBODY/
cp src/pair_chimes.{h,cpp} 		build/lammps_stable_29Oct2020/src/MANYBODY/
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
mv build/lammps_stable_29Oct2020/src/lmp_mpi_chimes exe

loc=`pwd`
echo ""
echo "Compilation complete. "
echo "Generated the following LAMMPS executable with ChIMES support:"
echo "${loc}/exe/lmp_mpi_chimes"
echo "See ${loc}/tests for usage examples"
echo ""



