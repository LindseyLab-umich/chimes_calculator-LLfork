set(CMAKE_Fortran_COMPILER "/opt/apps/gcc/13.2.0/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "13.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/opt/apps/gcc/13.2.0/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/opt/apps/gcc/13.2.0/bin/gcc-ar")
set(CMAKE_RANLIB "/opt/apps/gcc/13.2.0/bin/ranlib")
set(CMAKE_LINKER "/opt/apps/gcc/13.2.0/bin/ld")
set(CMAKE_Fortran_COMPILER_RANLIB "/opt/apps/gcc/13.2.0/bin/gcc-ranlib")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED TRUE)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/apps/gcc/13.2.0/lib/gcc/x86_64-pc-linux-gnu/13.2.0/finclude;/opt/intel/oneapi/ccl/2021.11/include;/opt/intel/oneapi/mpi/2021.11/include;/opt/intel/oneapi/mkl/2024.0/include;/opt/intel/oneapi/ipp/2021.10/include;/opt/intel/oneapi/ippcp/2021.9/include;/opt/intel/oneapi/dpl/2022.3/include;/opt/intel/oneapi/dpcpp-ct/2024.0/include;/scratch/projects/compilers/intel24.0/oneapi/dal/2024.0/include/dal;/opt/intel/oneapi/dev-utilities/2024.0/include;/opt/intel/oneapi/compiler/2024.0/include;/opt/intel/oneapi/tbb/2021.11/include;/opt/apps/gcc/13.2.0/lib/gcc/x86_64-pc-linux-gnu/13.2.0/include;/usr/local/include;/opt/apps/gcc/13.2.0/include;/opt/apps/gcc/13.2.0/lib/gcc/x86_64-pc-linux-gnu/13.2.0/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/apps/gcc/13.2.0/lib/gcc/x86_64-pc-linux-gnu/13.2.0;/opt/apps/gcc/13.2.0/lib64;/lib64;/usr/lib64;/opt/intel/oneapi/ccl/2021.11/lib;/opt/intel/oneapi/mpi/2021.11/opt/mpi/libfabric/lib;/opt/intel/oneapi/mpi/2021.11/lib;/opt/intel/oneapi/mkl/2024.0/lib;/opt/intel/oneapi/ipp/2021.10/lib;/opt/intel/oneapi/ippcp/2021.9/lib;/opt/intel/oneapi/dpl/2022.3/lib;/scratch/projects/compilers/intel24.0/oneapi/dal/2024.0/lib;/opt/intel/oneapi/compiler/2024.0/opt/compiler/lib;/opt/intel/oneapi/compiler/2024.0/lib;/opt/intel/oneapi/tbb/2021.11/lib;/opt/apps/gcc/13.2.0/x86_64-pc-linux-gnu/lib;/opt/apps/gcc/13.2.0/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
