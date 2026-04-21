#!/bin/bash
#SBATCH -J lmp_ins
#SBATCH -p skx-dev
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A TG-CHM250015
#SBATCH --output=build.out
#SBATCH --error=build.err

module purge
module load cmake
module load gcc
module load cuda/12.8
module load openmpi/5.0.8   

# Clean build
rm -rf build
mkdir -p build
cd build

# CRITICAL: Use nvcc_wrapper for proper CUDA+MPI compilation
NVCC_WRAPPER=$(realpath ../lib/kokkos/bin/nvcc_wrapper)
cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=on \
  -D BUILD_OMP=on \
  -D PKG_KOKKOS=on \
  -D Kokkos_ENABLE_CUDA=on \
  -D Kokkos_ENABLE_OPENMP=on \
  -D Kokkos_ARCH_HOPPER90=on \
  -D PKG_MANYBODY=on \
  -D PKG_MISC=on \
  -D PKG_SHOCK=on \
  -D PKG_REAXFF=on \
  -D PKG_SPIN=on \
  -D PKG_MOLECULE=on \
  -D PKG_CLASS2=on \
  -D PKG_KSPACE=on \
  -D CMAKE_CXX_COMPILER=$NVCC_WRAPPER \
  -D CMAKE_CXX_STANDARD=17

# Build
make -j 8

# Verify build
echo "=== Build complete. Checking executable ==="
ldd lmp | grep -E "(cuda|mpi)"