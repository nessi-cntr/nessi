cmake_minimum_required(VERSION 2.8)
project("libcntr")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# ~~ Commenly shared libs ~~
set(COMMON_LIBS "")

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

# ~~ Add HDF5 ~~
if (hdf5)
  #find_library(HDF5_LIB hdf5)
  find_package(HDF5)
  list(APPEND COMMON_LIBS ${HDF5_LIB})
  add_definitions("-DCNTR_USE_HDF5")
endif (hdf5)

# ~~ Add Eigen ~~
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

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
        set(MAIN_IN ${PROJECT_SOURCE_DIR}/doc/main.dox_inc)
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
