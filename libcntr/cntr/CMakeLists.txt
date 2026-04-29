# --  conditionally add ness source files to the library if ness=ON
set(NESS_FILES "")
if (ness)
    message(STATUS "Building with libness support")
    list(APPEND NESS_FILES ness.cpp ness2.cpp)
endif()


if (mpi) 
    # ~~ The actual target library ~~
    add_library(cntr SHARED
        fourier.cpp
        linalg_eigen.cpp
        integration.cpp
        integration_extern_templates.cpp
        cntr_convolution_extern_templates.cpp
        cntr_pseudo_convolution_extern_templates.cpp
        cntr_differentiation_extern_templates.cpp
        cntr_dyson_extern_templates.cpp
        cntr_dyson_omp_extern_templates.cpp
        cntr_equilibrium_extern_templates.cpp
        cntr_function_extern_templates.cpp
        cntr_herm_matrix_extern_templates.cpp
        cntr_herm_matrix_timestep_extern_templates.cpp
        cntr_herm_matrix_timestep_view_extern_templates.cpp
        cntr_herm_pseudo_extern_templates.cpp
        cntr_equilibrium_extern_templates.cpp
        cntr_utilities_extern_templates.cpp
        cntr_vie2_extern_templates.cpp
        cntr_distributed_array_extern_templates.cpp
        cntr_distributed_timestep_array_extern_templates.cpp
        cntr_getset_extern_templates.cpp
        cntr_mpitools_extern_templates.cpp
	# Conditionally add ness files if ness=ON
        ${NESS_FILES}
        )
else(mpi)
    # ~~ The actual target library ~~
    add_library(cntr SHARED
        fourier.cpp
        linalg_eigen.cpp
        integration.cpp
        integration_extern_templates.cpp
        cntr_convolution_extern_templates.cpp
        cntr_pseudo_convolution_extern_templates.cpp
        cntr_differentiation_extern_templates.cpp
        cntr_dyson_extern_templates.cpp
        cntr_dyson_omp_extern_templates.cpp
        cntr_equilibrium_extern_templates.cpp
        cntr_function_extern_templates.cpp
        cntr_herm_matrix_extern_templates.cpp
        cntr_herm_matrix_timestep_extern_templates.cpp
        cntr_herm_matrix_timestep_view_extern_templates.cpp
        cntr_herm_pseudo_extern_templates.cpp
        cntr_equilibrium_extern_templates.cpp
        cntr_utilities_extern_templates.cpp
        cntr_vie2_extern_templates.cpp
        cntr_getset_extern_templates.cpp
	# Conditionally add ness files if ness=ON
        ${NESS_FILES}
        )
    
endif(mpi)


set_property(TARGET cntr APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

# Link against FFTW if ness is enabled
		
# changed in NESSi 2.0 (libomp part)
if (hdf5 AND omp)
    add_subdirectory(hdf5)
    target_link_libraries(cntr cntr_hdf5 OpenMP::OpenMP_CXX)
elseif (hdf5 AND NOT omp)
    add_subdirectory(hdf5)
    target_link_libraries(cntr cntr_hdf5)
elseif (NOT hdf5 AND omp)
    target_link_libraries(cntr OpenMP::OpenMP_CXX)
endif()

# Always link in FFTW and other common libs
message(STATUS "linking cntr agaist ${COMMON_LIBS}")
target_link_libraries(cntr ${COMMON_LIBS})

# Add all headers as public-headers to `libcntr`. That way they will
# automatically be installed alongside `libcntr`.
file(GLOB cntr_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
set_target_properties(cntr PROPERTIES PUBLIC_HEADER "${cntr_HEADERS}")

# ~~ Installation targets ~~
install(TARGETS cntr
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/cntr)
if (hdf5)
    install(DIRECTORY hdf5 DESTINATION include/cntr
            FILES_MATCHING PATTERN "*.hpp")
endif (hdf5)
 install(DIRECTORY cntr/utils DESTINATION include/cntr
        FILES_MATCHING PATTERN "*.hpp")
install(FILES utils/read_inputfile.hpp DESTINATION include/cntr/utils)	


