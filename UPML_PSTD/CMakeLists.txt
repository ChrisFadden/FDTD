cmake_minimum_required (VERSION 2.8)

project(PSTD Fortran)
enable_language (Fortran)
set(MY_SOURCES src/Parameters.f90 src/Grid.f90 
src/UPML.f90 src/FFTW_utils.f90 src/MPI_utils.f90 src/Update.f90 src/PrintField.f90 
 src/main.f90)
           
set(CMAKE_Fortran_COMPILER /usr/bin/mpifort)

find_library(fftw fftw3)
find_path(fftw_path fftw3.h)

find_package(MPI REQUIRED)
find_package(OpenMP)
link_directories(${MPI_Fortran_LIBRARIES})
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_EXE_LINKER_FLAGS "-lfftw3 -lpthread -fopenmp")
set(CMAKE_FORTRAN_FLAGS_DEBUG "-g -Wall -fbounds-check ffpe-trap=invalid,zero,overflow")

set(CMAKE_FORTRAN_COMPILE_FLAGS ${CMAKE_FORTRAN_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ${OpenMP_C_FLAGS})
include_directories(${fftw_path} ${MPI_INCLUDE_PATH})

add_executable(PSTD ${MY_SOURCES})
set(CMAKE_BUILD_TYPE Debug)
set(ENABLE_MPI On)
