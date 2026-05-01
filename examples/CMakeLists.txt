cmake_minimum_required(VERSION 3.10)
project("demo_programs")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(COMMON_LIBS "")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/exe)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/exe)

# Optionally enable OpenMPI, OpenMP, HDF5

option(ness "Enable NESS support" OFF)
#if (ness)
#   message(STATUS "Building with NESS support")

    #  Find FFTW
    #  find_library(FFTW3_LIB fftw3 PATHS ${CMAKE_LIBRARY_PATH} NO_DEFAULT_PATH)
    #  find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h PATHS ${CMAKE_INCLUDE_PATH})

    #if (FFTW3_LIB AND FFTW3_INCLUDE_DIR)
      #  message(STATUS "FFTW library found: ${FFTW3_LIB}")
       # message(STATUS "FFTW include directory found: ${FFTW3_INCLUDE_DIR}")
        #include_directories(${FFTW3_INCLUDE_DIR})
       # add_definitions("-DCNTR_USE_FFTW")
       # list(APPEND COMMON_LIBS ${FFTW3_LIB})
    #else()
      #  message(FATAL_ERROR "FFTW library or headers not found")
    #endif()
#endif()

option(hdf5 "Build with HDF5 support" OFF)
if (hdf5)
   message(STATUS "Building with HDF5")
   find_library(HDF5_LIB hdf5)
   if (HDF5_LIB)
      message(STATUS "Found HDF5: ${HDF5_LIB}")
      list(APPEND COMMON_LIBS ${HDF5_LIB})
   else()
      message(FATAL_ERROR "HDF5 library not found!")
   endif()
   add_definitions("-DCNTR_USE_HDF5")
endif (hdf5)

# ~~ OpenMP and MPI flags ~~

option(omp "Build with OpenMP support" OFF)
if (omp)
    message(STATUS "Building with OpenMP")
    find_package(OpenMP REQUIRED)
    add_definitions("-DCNTR_USE_OMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

option(mpi "Build with OpenMPI support" OFF)

if (mpi)
    message(STATUS "Building with OpenMPI")
    find_package(MPI REQUIRED)
    add_definitions("-DCNTR_USE_MPI")
    set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_COMPILE_FLAGS
        "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS
        "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND COMMON_LIBS ${MPI_LIBRARIES})
endif ()

# ~~ Add Eigen ~~
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

add_subdirectory(programs)


# ~~ Include Guard ~~
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")
