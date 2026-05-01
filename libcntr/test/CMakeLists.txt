# -- Catch tests compiled into target "runtest"

# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Base test sources
set(TEST_SOURCES
    runtest.cpp
    bubble.cpp    
    convolution.cpp
    downfold.cpp
    dyson.cpp
    dyson_new.cpp
    dyson_moving.cpp
    moving_data_exchange.cpp
    moving_green_textio.cpp
    moving_green_simple_algebra.cpp	               
    equilibrium.cpp
    function.cpp
    getset.cpp
    getset_nonherm.cpp
    green_single_pole_XX.cpp
    herm_matrix_algebra.cpp
    herm_matrix_member.cpp  
    herm_matrix_readwrite.cpp
    herm_matrix_setget_timestep.cpp
    herm_matrix_submatrix.cpp
    herm_member_timestep.cpp
    herm_member_timestep_view.cpp
    integration.cpp
    linalg.cpp
    matsubara.cpp    
    utilities.cpp
)

# Add tests using hdf5 if ON
if (hdf5)
    message(STATUS "Including hdf5 test files")

    list(APPEND TEST_SOURCES
    moving_green_hdf5.cpp
    fkm_bethe_quench.cpp
    herm_matrix_hdf5.cpp
    herm_matrix_readwrite_hdf5.cpp
    )
endif()

# Add Ness2-specific tests if ness is ON
if (ness)
    message(STATUS "Including Ness2 test files")

    list(APPEND TEST_SOURCES
        ness2_fft_array_basic.cpp
        ness2_fft_array_elements.cpp
        ness2_fft_array_ops.cpp
        ness2_herm_matrix_ness_basic.cpp
        ness2_herm_matrix_ness_elements.cpp
        ness2_herm_matrix_ness_ops.cpp
	ness2_dyson.cpp
	ness2_grid_info.cpp
	ness2_omp_manager.cpp
	ness2_herm_matrix_ness_integral_transform.cpp
	ness2_resample.cpp
	ness2_density_matrix.cpp
    )
endif()

# Build runtest executable
add_executable(runtest ${TEST_SOURCES})

# Link libraries
target_link_libraries(runtest cntr ${COMMON_LIBS})

# MPI variant
if (mpi)
    add_executable(runtest_mpi runtest_mpi.cpp)
    target_link_libraries(runtest_mpi cntr ${COMMON_LIBS})
endif()

# -- "test_all" target for running all the test using Catch directly
add_custom_target(test COMMAND runtest -d yes DEPENDS runtest)

if (mpi)
    add_custom_target(test_mpi COMMAND ${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/runtest_mpi -d yes DEPENDS runtest_mpi)
endif()

