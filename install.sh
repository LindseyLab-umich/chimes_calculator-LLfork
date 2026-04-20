#!/bin/bash

# Builds all relevant chimes_calculator executables/library files 
# Run with:
# ./install.sh 
# or
# ./install.sh <debug option (0 or 1)> <install prefix (full path)>  <LLNL computer (0 or 1)> <VERBOSE 0 or 1>

# Verbose is required for testing.

DEBUG=${1-0}  # False (0) by default.
PREFX=${2-""} # Empty by default
LLNL=${3-0}
VERBOSE=${4-0}

echo "DEBUG=$DEBUG"
echo "PREFIX=$PREFX"
echo "LLNL=$LLNL"
echo "VERBOSE=$VERBOSE"

# Clean up previous installation, 

./uninstall.sh $PREFX

# Move into build directory 

mkdir build
cd build

if [[ $VERBOSE -eq 1 ]] ; then
    my_flags="-DVERBOSE=1"
else
    my_flags="-DVERBOSE=0"
fi    

# Generate cmake flags
#
# CMAKE is broken for fortran builds with Intel compilers.
# Uncomment to use intel compilers for C and C++ interfaces.
#
if [[ $LLNL -eq 1 ]] ; then
    if [[ "$SYS_TYPE" == "toss_3_x86_64_ib" ]] ; then
	# Try newest intel compiler for best performance.
	module load cmake/3.14.5
	module load intel/2021.3
	ICC=`which icc`
	IFORT=`which ifort`
	my_flags="$my_flags -DCMAKE_CXX_COMPILER=${ICC} -DCMAKE_Fortran_COMPILER=${IFORT} -DCMAKE_C_COMPILER=${ICC}"
	my_flags="$my_flags -DCMAKE_CXX_FLAGS_RELEASE=\"-O3 -fno-alias -fno-fnalias -xhost\""
    elif [[ "$SYS_TYPE" == "toss_4_x86_64_ib" ]] ; then
	module load cmake/3.14.5
	module load intel-classic/19.1.2
	ICC=`which icc`
	IFORT=`which ifort`
	my_flags="$my_flags -DCMAKE_CXX_COMPILER=${ICC} -DCMAKE_Fortran_COMPILER=${IFORT} -DCMAKE_C_COMPILER=${ICC}"
	my_flags="$my_flags -DCMAKE_CXX_FLAGS_RELEASE=\"-O3 -fno-alias -fno-fnalias -xhost\""
    else	
	echo "Unknown LLNL operating system type: $SYS_TYPE"
    fi
else
    echo "Compiling for a non-LLNL computer system"
fi    

if [ ! -z $PREFX ] ; then
	my_flags="-DCMAKE_INSTALL_PREFIX=${PREFX} $my_flags"
fi

if [ $DEBUG -eq 1 ] ;then

	my_flags="${my_flags} -DDEBUG=1 -DCMAKE_BUILD_TYPE=Release" 
else
	my_flags="${my_flags} -DDEBUG=0 -DCMAKE_BUILD_TYPE=Release" 
fi

# Setup, make and install

echo "Invoking cmake $my_flags"
cmake $my_flags ..
make
if [ ! -z $PREFX ] ; then
	make install
fi

cd ..
