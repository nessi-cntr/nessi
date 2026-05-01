cmake_minimum_required(VERSION 3.10)
project("libcntr")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


# ~~ Only relevant for macOS ~~
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()
set(CMAKE_MACOSX_RPATH 1)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building in release mode")
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -DNDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")


# Set the runtime path where the executable will look for libraries
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

# Optionally: Set install name directory if you want specific control (less commonly needed)
# set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)


# ~~ Commenly shared libs ~~
set(COMMON_LIBS "")

# ~~ OpenMP and MPI flags ~~
option(omp "Build with OpenMP support" OFF)
if (omp)
  message(STATUS "Building with OpenMP")
  find_package(OpenMP REQUIRED)

  if(OpenMP_CXX_FOUND)
    list(APPEND COMMON_LIBS OpenMP::OpenMP_CXX)
  endif()
  if(OpenMP_C_FOUND)
    list(APPEND COMMON_LIBS OpenMP::OpenMP_C)
  endif()

  add_definitions("-DCNTR_USE_OMP")

  # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (omp)



option(mpi "Build with OpenMPI support" OFF)
if (mpi)
    message(STATUS "Building with OpenMPI")
    find_package(MPI REQUIRED)

    message(STATUS "MPI_C_COMPILER: ${MPI_C_COMPILER}")
    message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
    message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
    message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
    message(STATUS "MPI_COMPILE_FLAGS: ${MPI_COMPILE_FLAGS}")
    message(STATUS "MPI_LINK_FLAGS: ${MPI_LINK_FLAGS}")

    add_definitions("-DCNTR_USE_MPI")
    set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND COMMON_LIBS ${MPI_LIBRARIES})
endif (mpi)

# ~~ Add HDF5 ~~
option(hdf5 "Build with HDF5 support" OFF)
if (hdf5)
  message(STATUS "Building with HDF5 support")

  if(hdf5_prefix)
    find_library(HDF5_LIB
        NAMES hdf5
        HINTS "${hdf5_prefix}/lib"
    )    
    include_directories(${hdf5_prefix}/include)
    link_directories(${hdf5_prefix}/lib)    
  else(hdf5_prefix)
    find_package(HDF5)
    set(HDF5_LIB ${LIBS} ${HDF5_LIBRARIES})
    include_directories(${HDF5_INCLUDE_DIR})
  endif(hdf5_prefix)

  list(APPEND COMMON_LIBS ${HDF5_LIB})
  add_definitions("-DCNTR_USE_HDF5")
endif (hdf5)

# ~~ Add Eigen ~~
find_package(Eigen3 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    message(STATUS "Eigen3 successfully found!")
    message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "Eigen3 not found! Please install Eigen3 or set EIGEN3_INCLUDE_DIR manually.")
endif()
list(APPEND COMMON_LIBS Eigen3::Eigen)
#include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# ~~ Add FFTW (required for NESS) ~~
option(ness "Build with NESS support" OFF)
if (ness)
    message(STATUS "NESS enabled")
    add_definitions(-DCNTR_USE_NESS)

    if(fftw3_prefix) 
        find_library(FFTW_LIB
            NAMES fftw3
            HINTS "${fftw3_prefix}/lib"
        )     

        link_directories(${fftw3_prefix}/lib)
        include_directories(${fftw3_prefix}/include)

        # Try OpenMP version first
        if (omp) 
            find_library(FFTW_OMP_LIB
                NAMES fftw3_omp
                HINTS "${fftw3_prefix}/lib"
            )        

            if (FFTW_OMP_LIB)
                message(STATUS "FFTW-OMP library found: ${FFTW_OMP_LIB}")
                list(APPEND COMMON_LIBS ${FFTW_OMP_LIB})
                add_definitions(-DNESS2_FFT_USE_OMP=1)         # directory-scope macro
            else (FFTW_OMP_LIB)

                # Try pthreads version next
                find_library(FFTW_THREADS_LIB
                    NAMES fftw3_threads
                    HINTS "${fftw3_prefix}/lib"
                )        

                find_package(Threads)
                message(STATUS "Thread libs: ${CMAKE_THREAD_LIBS_INIT}")
    
                if (FFTW_THREADS_LIB AND Threads_FOUND)
                    message(STATUS "FFTW pthreads library found: ${FFTW_THREADS_LIB}")
                    list(APPEND COMMON_LIBS ${FFTW_THREADS_LIB} Threads::Threads)
                    add_definitions(-DNESS2_FFT_USE_OMP=1)         # directory-scope macro
                else ()
                    message(WARNING "omp=ON but neither fftw3_omp nor fftw3_threads found; using single-threaded FFTW")

                    find_library(FFTW_LIB
                            NAMES fftw3
                            HINTS "${fftw3_prefix}/lib"
                        ) 
                        
                    if (FFTW_LIB)
                        message(STATUS "Using single-threaded FFTW library: ${FFTW_LIB}")
                        list(APPEND COMMON_LIBS ${FFTW_LIB})
                    else (FFTW_LIB)
                        message(FATAL_ERROR "FFTW library not found in ${fftw3_prefix}/lib")
                    endif (FFTW_LIB)

                endif ()
            endif (FFTW_OMP_LIB)

        else (omp)
            message(STATUS "omp=OFF: using single-threaded FFTW")

            find_library(FFTW_LIB
                        NAMES fftw3
                        HINTS "${fftw3_prefix}/lib"
                    )        

            if (FFTW_LIB)
                message(STATUS "Using single-threaded FFTW library: ${FFTW_LIB}")
                list(APPEND COMMON_LIBS ${FFTW_LIB})
            else (FFTW_LIB)
                message(FATAL_ERROR "FFTW library not found in ${fftw3_prefix}/lib")
            endif (FFTW_LIB)

        endif (omp)

    else(fftw3_prefix)
        message(STATUS "Searching for FFTW3 in system paths")

        find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)
        find_library(FFTW3_LIB fftw3 REQUIRED)

        if (NOT FFTW3_LIB)
            message(FATAL_ERROR "FFTW3 library not found")
        endif()

        include_directories(${FFTW3_INCLUDE_DIR})
        list(APPEND COMMON_LIBS ${FFTW3_LIB})

        if (omp) 
            # try to find fftw3_omp first
            find_library(FFTW_OMP_LIB fftw3_omp)

            if(FFTW_OMP_LIB)
                message(STATUS "FFTW-OMP library found: ${FFTW_OMP_LIB}")
                list(APPEND COMMON_LIBS ${FFTW_OMP_LIB})
                add_definitions(-DNESS2_FFT_USE_OMP=1)         # directory-scope macro
            else (FFTW_OMP_LIB)
                # try to find fftw3_threads next
                find_library(FFTW3_THREADS_LIB fftw3_threads)
                find_package(Threads)
                message(STATUS "Thread libs: ${CMAKE_THREAD_LIBS_INIT}")

                if(FFTW3_THREADS_LIB AND Threads_FOUND)
                    message(STATUS "Using FFTW pthreads backend: ${FFTW3_THREADS_LIB}")
                    list(APPEND COMMON_LIBS ${FFTW3_THREADS_LIB} Threads::Threads)
                    add_definitions(-DNESS2_FFT_USE_OMP=1)         # directory-scope macro
                else()
                    message(WARNING "omp=ON but neither fftw3_omp nor fftw3_threads found; using single-threaded FFTW")

                    find_library(FFTW3_LIB fftw3 REQUIRED)
                    if (FFTW3_LIB)
                        message(STATUS "Using single-threaded FFTW library: ${FFTW3_LIB}")
                        list(APPEND COMMON_LIBS ${FFTW3_LIB})
                    else (FFTW3_LIB)
                        message(FATAL_ERROR "FFTW library not found")
                    endif(FFTW3_LIB)

                endif()
            endif (FFTW_OMP_LIB)

        else (omp)
            message(STATUS "omp=OFF: using single-threaded FFTW")

            find_library(FFTW3_LIB fftw3 REQUIRED)
            if (FFTW3_LIB)
                message(STATUS "Using single-threaded FFTW library: ${FFTW3_LIB}")
                list(APPEND COMMON_LIBS ${FFTW3_LIB})
            else (FFTW3_LIB)
                message(FATAL_ERROR "FFTW library not found")
            endif(FFTW3_LIB)

        endif (omp)
    endif (fftw3_prefix)

endif(ness)


# ~~ Relocatable code for shared library ~~
set(GCC_RELOCATABLE_FLAGS "-fPIC")
add_definitions(${GCC_RELOCATABLE_FLAGS})

# ~~ Warning Flags ~~
# Re-enable whenever you feel strong enough.
# set(GCC_WARNING_FLAGS "-Wall -Wextra -pedantic")
# add_definitions(${GCC_WARNING_FLAGS})

# ~~ CMake Configuration ~~
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cntr/)

add_subdirectory(cntr)

option(BUILD_TEST "Build test suite" ON)
if(BUILD_TEST)
    add_subdirectory(test)
endif(BUILD_TEST)

# ~~ Generation of documentation ~~
# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

if(BUILD_DOC)
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/Doxyfile.txt)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        #set(DOXYGEN_OUT ${PROJECT_SOURCE_DIR}/Doxyfile.txt)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message(STATUS "Doxygen build started")

         # set input and output files
        set(MAIN_IN ${PROJECT_SOURCE_DIR}/main.dox_inc)
        set(MAIN_OUT ${CMAKE_CURRENT_BINARY_DIR}/main.dox)
        configure_file(${MAIN_IN} ${MAIN_OUT} @ONLY)

        # copy images
        file(GLOB IMAGES ${CMAKE_CURRENT_SOURCE_DIR}/doc/*.png)
        file(COPY ${IMAGES}
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/html/)
        

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
      message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif(BUILD_DOC)

# ~~ Include Guard ~~
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")
