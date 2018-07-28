#src/CMakeLists.txt

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

set(EXE Main)
MESSAGE( "\n The name of the executable is: ${EXE} \n")
# FIND LIBRARIES

find_library(MPI_LIBRARY NAMES mpi )
MESSAGE( "mpi library variable is: ${MPI_LIBRARY}")

#find_library(BLAS_LIBRARY NAMES libblas HINTS "/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/libblas.so")
#set (BLAS_LIBRARY "/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/libblas.so")

#find_library(Blas NAMES blas HINTS "/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/libblas.so")
#MESSAGE( "blas library variable is: ${Blas}")


find_library(LAPACK_LIBRARY NAMES lapack HINTS "/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/")
MESSAGE( "lapack library variable is: ${LAPACK_lIBRARY}")

find_library(PARMETIS_LIBRARY NAMES parmetis HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../lib/")
MESSAGE( "parmetis library variable is: ${PARMETIS_LIBRARY}")

find_library(HDF5_LIBRARY NAMES hdf5 HINTS "/u/sw/")
MESSAGE( "hdf5 library variable is: ${HDF5_LIBRARY}")

find_library(QHULL_LIBRARY NAMES qhull HINTS "/u/sw/")
MESSAGE( "qhull library variable is: ${QHULL_LIBRARY}")

find_library(TRILINOS_LIBRARY NAMES trilinos epetra)
MESSAGE( "trilinos library variable is: ${TRILINOS_LIBRARY}")

#FIND_PACKAGE(TRILINOS REQUIRED  )
#MESSAGE(${TRILINOS_FIND_VERSION})

find_library(MUPARSER_LIBRARY NAMES muparser HINTS "/u/sw/")
MESSAGE( "muparser library variable is: ${MUPARSER_LIBRARY}")

#find_library(TEUCHOS_LIBRARY NAMES teuchos HINTS "/u/sw/")
#MESSAGE( "teuchos library variable is: ${TEUCHOS_LIBRARY}")

INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/include/)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/lib/)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/lifev)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/lifev/static_lib)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/parmetis/include)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/parmetis/lib)





#set(Openblas_DIR /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/lib/cmake)
#message("\n\nopenblas dir: ${Openblas_DIR}")
find_package(Openblas 0.2.17 NAMES OpenBLAS CONFIG)
if(Openblas_FOUND)
  MESSAGE( "Openblas found:")
  MESSAGE( "              : ${Openblas_INCLUDE_DIRS}")
  MESSAGE( "              : ${Openblas_LIBRARY_DIRS}")
  include_directories(${Openblas_INCLUDE_DIRS})
  link_directories(${Openblas_LIBRARY_DIRS})
endif()
message("\n\nopenblas dir: ${Openblas_DIR}\n\n")



find_package(Lapack 3.6.0 NAMES lapack  HINTS /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/cmake)
if(Lapack_FOUND)
  MESSAGE( "Lapack found:")
  MESSAGE( "              : ${Lapack_INCLUDE_DIRS}")
  MESSAGE( "              : ${Lapack_LIBRARY_DIRS}")
  include_directories(${Lapack_INCLUDE_DIRS})
  link_directories(${Lapack_LIBRARY_DIRS})
endif()

find_package(Trilinos 12.6.3)
if(Trilinos_FOUND)
  MESSAGE( "Trilinos found:")
  MESSAGE( "              : ${Trilinos_INCLUDE_DIRS}")
  MESSAGE( "              : ${Trilinos_LIBRARY_DIRS}")
  include_directories(${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
endif()

find_package(Boost 1.60.0)
if(Boost_FOUND)
  MESSAGE( "BOOST found ${Boost_LIBRARIES}")
  include_directories(${Boost_INCLUDE_DIRS})
endif()

set(BLAS_DYN_LIB /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64/libblas.so)

set(BELOSEPETRA_DYN_LIB /u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/trilinos/12.6.3/lib/libbelosepetra.so)



INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/muparser/2.2.5/include)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/muparser/2.2.5/lib)



INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/muparser/2.2.5/include)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/muparser/2.2.5/lib)

INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hdf5/1.8.16/include)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hdf5/1.8.16/lib)

INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/eigen/3.2.8/include/eigen3)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/eigen/3.2.8/lib)

INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hypre/2.11.0/include)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/hypre/2.11.0/lib)

INCLUDE_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/include)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/lib)

#LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/lib)
#LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/base/lib)

LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/lapack/3.6.0/lib64)
#LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/metis/5/lib)


set (corelib         "${CMAKE_SOURCE_DIR}/lib/lifev/static_lib/liblifevcore.a")
set (etalib          "${CMAKE_SOURCE_DIR}/lib/lifev/static_lib/liblifeveta.a")
set (navierstokeslib "${CMAKE_SOURCE_DIR}/lib/lifev/static_lib/liblifevnavierstokes.a")
set (himodlib        "${CMAKE_SOURCE_DIR}/lib/lifev/static_lib/liblifevhimod.a")

#set(LIFEV_PATH "/home/zampieri/Desktop/FSI_HiMod/cmake_project/lib/lifev/" )
#find_package( LifeV PATHS ${LIFEV_PATH}/lib/cmake/LifeV ${LIFEV_PATH} )
#MESSAGE(${LIFEV_FIND_VERSION})

add_definitions(-std=c++11)
SET(GCC_COVERAGE_LINK_FLAGS)
set(CMAKE_CXX_STANDARD 11)
add_definitions(${GCC_COVERAGE_LINK_FLAGS})

file(GLOB sourcefiles "*.cpp")
add_executable(${EXE} ${sourcefiles} ${headerfiles})


#file(GLOB sourcefiles "test.cpp")
#add_executable(${EXE} ${sourcefiles})

target_link_libraries(${EXE}
        #${TRILINOS_LIBRARY}
        blas muparser
        belosepetra
        teuchosnumerics
        belos
        epetra
        lifevcore lifeveta lifevnavierstokes lifevhimod
        locathyra locaepetra  localapack  loca  noxepetra  noxlapack  nox
         rythmos  teko  ifpack2-adapters  ifpack2  amesos2
         stratimikos  stratimikosbelos  stratimikosaztecoo  stratimikosamesos  stratimikosml
         stratimikosifpack  shylu  ml  galeri-xpetra
         galeri-epetra  isorropia  xpetra-sup  xpetra  thyratpetra  thyraepetraext
         ifpack  amesos  anasazitpetra  ModeLaplace  anasaziepetra
         anasazi  thyraepetra  thyracore  rtop  belostpetra  belosepetra  belos
         tpetraext  tpetrainout  tpetra  kokkostsqr  tpetrakernels
         tpetraclassiclinalg  tpetraclassicnodeapi  tpetraclassic  kokkosalgorithms  kokkoscontainers
         aztecoo  zoltan  epetraext  triutils  epetra
         teuchoskokkoscomm  teuchoskokkoscompat  teuchosremainder  teuchosnumerics  teuchoscomm
         teuchosparameterlist  teuchoscore  kokkoscore
        ${PARMETIS_LIBRARY}
        ${HDF5_LIBRARY}
        #${BLAS_LIBRARY}
        ${MPI_LIBRARY}
        ${BLAS_DYN_LIB}
        ${BELOSEPETRA_DYN_LIB}
        ${corelib}
        ${etalib}
        ${navierstokeslib}
        ${himodlib}
        #${TRILINOS_LIBRARY}

        #
        )

MESSAGE("\n Finished the source CMakeLists \n")
